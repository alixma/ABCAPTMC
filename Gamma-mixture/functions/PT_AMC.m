function [X, n, nsw] = PT_AMC(K, T, N, C, lambdas, Lambda, rho, alpha, theta, p, init, npairs, deltat, correct)
% Single processor Anytime Parallel Tempeting Monte Carlo algorithm
% Note: Only outputs the cold chain
% Inputs:   K: # chains
%           T: units of virtual time for which to run algorithm
%           N:  # slots allocated to matrix X to record samples
%           C:  run the algorithm C times in parallel
%           lambdas, Lambda: vector of temperatures, # temperatures (usually Lambda = W)
%           rho: standard deviation for Metropolis update
%           alpha, theta: parameters of Gamma mixture
%           p: computational complexity
%           init: point of initialisation of chains
%           npairs # pairs of adjacent chains to swap in exchange moves
%           deltat: # time steps between exchange moves
%           correct: apply bias correction (1) or not (0)
% Outputs:  X: cold chain(s)
%           n: sample size obtained for each chain 
%           nsw: total # exchange moves performed   

if ~exist('init', 'var')
    init='prior';
end

tis = (1+deltat):deltat:T; 
X = zeros(N, C);
n = zeros(K, C);
nsw = 0;

if (correct)
    even=1:2:(K-2);
    odd=2:2:(K-2);
    pair=[1:(K-2);  2:(K-1)]';
else
    even=1:2:(K-1);
    odd=2:2:(K-1);
    pair=[1:(K-1);  2:K]';
end
% tic
parfor c=1:C
    X1 = zeros(N, 1);
    evenodd=0;
    pair1=pair;
    lbds=lambdas;
    
    % initialise K chains
    if strcmp(init, 'prior')
        x = gamrnd_mixture(alpha, theta, K, 1);
    else
        x = zeros(K, 1) + init;
    end
    
    % simulate chains one by one
    k = 1;
    ti = 1;
    update_first=0;

    h = zeros(K,1);
    n1 = ones(K,1);
    X1(1) = x(1);

    for t = tis
        %% LOCAL MOVES %%
        % simulate chains one by one
        while ti <= t
            
            %   update k-th state directly if previous hold time ended  
            %at the same time as the deadline
            if(~update_first) 
                % hold k-th state
                if(h(k)<=1)
                    h(k) = h(k) + gamrnd_mixture(x(k)^p./theta, theta, 1, 1, 1);
                end
                while (h(k) > 1) && (ti <= t)
                    ti = ti + 1;
                    h(k) = h(k) - 1;
                end
                if(h(k) > 1)
                    update_first=1;
                end
            end
            if (ti <= t)  %do not update if deadline has occured 
                n1(k) = n1(k) + 1;
                % update k-th state using random walk Metropolis
                can = normrnd(x(k), rho);
                % odds ratio
                A1 = gampdf_mixture_temp(can, lbds(k), Lambda, alpha, theta, 0.5);
                A2 = gampdf_mixture_temp(x(k), lbds(k), Lambda, alpha, theta, 0.5);
                A = A1/A2;
                % accept proposal with probability A
                if(unifrnd(0,1) < A)
                    x(k) = can;
                end
                if(k==1)
                    X1(n1(k)) = x(k);
                end
                update_first=0;
                
                % next chain
                k = k + 1;
                if (k > K)
                    k = 1;
                end
            end
        end
        
        %% EXCHANGE MOVES %%
        % perform exchange moves before resuming local moves
        nsw=nsw+1;
        if(correct)
            % get rid of k-th chain
            exK =  [1:(k-1) (k+1):K];
        else
            exK = 1:K;
        end
        XE = x(exK);
        nE = n1(exK);
        lbdsE = lbds(exK);
        if(~correct)
            hE = h(exK); 
        end
        
        % perform exchange moves on remaining chains
        if(npairs < min(length(even),length(odd)))
            % either select a subset of even/odd pairs
            if(evenodd)
                swaps=randsample(odd, npairs, false);
                evenodd = 0;
            else
                swaps=randsample(even, npairs, false);
                evenodd = 1;
            end
        else % or select all even/odd pairs (default)
            if(evenodd)
                swaps=odd;
                evenodd = 0;
            else
                swaps=even;
                evenodd = 1;
            end
        end
        nE(swaps)=nE(swaps)+1;
        nE(swaps+1)=nE(swaps+1)+1;
        n1(exK) = nE;
        % swapping all selected chains at once
        swap = pair1(swaps,:);
        i=swap(:,1);
        j=swap(:,2);
        % odds ratio
        A1=gampdf_mixture_temp(XE(i)', lbdsE(j), Lambda, alpha, theta).*gampdf_mixture_temp(XE(j)', lbdsE(i), Lambda, alpha, theta);
        A2=gampdf_mixture_temp(XE(i)', lbdsE(i), Lambda, alpha, theta).*gampdf_mixture_temp(XE(j)', lbdsE(j), Lambda, alpha, theta);
        A=A1./A2;
        % accept or reject swaps
        mkswap = unifrnd(0,1,size(A)) < A;
        % perform the accepted swaps
        x2 = XE(swap(mkswap,2));
        XE(swap(mkswap,2)) = XE(swap(mkswap,1));
        XE(swap(mkswap,1)) = x2;        
        if(~correct) % swap hold times too if no bias correction
            h2 = hE(swap(mkswap,2));
            hE(swap(mkswap,2)) = hE(swap(mkswap,1));
            hE(swap(mkswap,1)) = h2;
            h(exK) = hE;
        end
        
        %record exchanges on chains
        x(exK) = XE;
        X1(n1(1)) = x(1);
    end
    X(:,c) = X1;
    n(:,c) = n1;
end
% toc
maxn = max(n(1,:));
X = X(1:maxn,:);
end
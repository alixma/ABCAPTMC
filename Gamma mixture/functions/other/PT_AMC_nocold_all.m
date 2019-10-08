function [X, Ac, Sw, n, nsw] = PT_AMC_nocold_all(K, T, N, C, lambdas, Lambda, rho, alpha, theta, p, init, npairs, deltat, correct)
% Single processor Anytime Parallel Tempeting Monte Carlo algorithm
% Notes:    No local moves on the cold chain
%           Outputs all chains
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
% Outputs:  X: all chains
%           Ac: accepted local moves for each chain
%           Sw: accepted swaps for each chain
%           n: sample size obtained for each chain
%           nsw: total # exchange moves performed

if ~exist('init', 'var')
    init='prior';
end

tis = (1+deltat):deltat:T; 
X = zeros(K, N, C);
Ac = zeros(K, N, C);
Sw = zeros(K, size(tis,2), C);
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
%tic
parfor c=1:C
    X1 = zeros(K, N);
    Ac1 = zeros(K, N);
    Sw1 = zeros(K, size(tis,2));   
    pair1=pair
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
    X1(:, 1) = x;
         
    for t = tis
        %% LOCAL MOVES %%
        % simulate chains one by one
        while ti <= t
            if(k ~= 1) % do not hold the cold chain
                %   update k-th state directly if previous hold time ended
                %at the same time as the deadline
                if(~update_first)
                    % hold k-th state
                    if(h(k)<=1)
                        h(k) = h(k) + gamrnd_mixture(x(k)^p./theta, theta, 1, 1, 1);
                    end
                    while h(k) > 1 && ti <= t
                        ti = ti + 1;
                        h(k) = h(k) - 1;
                    end
                    if(h(k) > 1)
                        update_first=1;
                    end
                end
            end
            if (ti <= t) %do not update if deadline has occured
                if(k ~= 1) % do not update the cold chain
                    n1(k) = n1(k) + 1;
                    % update k-th state using random walk Metropolis
                    can = normrnd(x(k), rho);
                    % odds ratio
                    A1 = gampdf_mixture_temp(can, lbds(k), Lambda, alpha, theta);
                    A2 = gampdf_mixture_temp(x(k), lbds(k), Lambda, alpha, theta);
                    A = A1/A2;
                    % accept proposal with probability A
                    if(unifrnd(0,1) < A)
                        x(k) = can;
                        Ac1(k, n1(k))=1;
                    end
                    X1(k, n1(k)) = x(k);
                    update_first=0;
                end
                
                % next chain
                k = k + 1;
                if (k > K)
                    k = 1;
                end
            end
        end
                 
        %% EXCHANGE MOVES %%
        % perform exchange moves before resuming update
        nsw=nsw+1;
        if(correct)
            % get rid of k-th chain
            exK =  [1:(k-1) (k+1):K];
            L = K-1;
        else
            exK = 1:K;
            L = K;
        end
        XE = zeros(L,1);
        SwE = zeros(L,1);
        AcE = zeros(L,1);
        if(~correct)
            hE = h(exK); 
        end
        lbdsE=lbds(exK);
        
        for ll = 1:L
            XE(ll) = X1(exK(ll), n1(exK(ll)));
        end
        
        % perform exchange moves on remaining chains
        if(npairs < min(length(even),length(odd)))
            % either select a subset of even/odd pairs
            if(mod((T-t)./deltat, 2)==0)
                swaps=randsample(odd, npairs, false);
            else
                swaps=randsample(even, npairs, false);
            end
        else % or select all even/odd pairs (default)
            if(mod((T-t)./deltat, 2)==0)
                swaps=odd;
            else
                swaps=even;
            end
        end
        n1(exK(swaps))=n1(exK(swaps))+1;
        n1(exK(swaps+1))=n1(exK(swaps+1))+1;
        
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
        if(~correct) %swap hold times too if no bias correction
            h2 = hE(swap(mkswap,2));
            hE(swap(mkswap,2)) = hE(swap(mkswap,1));
            hE(swap(mkswap,1)) = h2;
            h(exK) = hE;
        end

        % record accepted swaps
        SwE(swap(mkswap,:)) = 1;
        AcE(swap(mkswap,:)) = 2;
        
        %record exchanges on chains
        for ll = 1:L
            X1(exK(ll), n1(exK(ll))) = XE(ll);
            x(exK(ll))= XE(ll);
            Sw1(exK(ll), (ti-2)/deltat) =  SwE(ll);
            Ac1(exK(ll), n1(exK(ll))) = AcE(ll);
        end
    end
    X(:,:,c) = X1;
    Ac(:,:,c) = Ac1;
    Sw(:,:,c) = Sw1;
    n(:,c) = n1;
end
%toc
maxn = max(n(:));
X = X(:,1:maxn,:);
Ac = Ac(:,1:maxn,:);
end
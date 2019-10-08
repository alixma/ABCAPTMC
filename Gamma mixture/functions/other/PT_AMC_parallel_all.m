function [Xout, n] = PT_AMC_parallel_all(W, Kk, T, N, lambdas, Lambda, rho, alpha, theta, p, init, npairs, deltat, correct)
% Multi-processor  Anytime Parallel Tempering Monte Carlo algorithm
% Notes:    Outputs all chains
% Inputs:   W:  # workers
%           Kk: # chains per worker
%           T: units of virtual time for which to run algorithm
%           N:  # slots allocated to matrix X to record samples
%           lambdas, Lambda: vector of temperatures, # temperatures (usually Lambda = W)
%           rho: standard deviation for Metropolis update
%           alpha, theta: parameters of Gamma mixture
%           p: computational complexity
%           npairs: # pairs of adjacent chains to swap in exchange moves
%           deltat: # time steps between exchange moves
%           correct: apply bias correction (1) or not (0)
% Outputs   Xout: all chains
%           n: sample size obtained for each chain

if ~exist('init', 'var')
    init='prior';
end

Xout = zeros(Kk, N, Lambda);
tis = (1+deltat):deltat:T; % real-time schedule
n = ones(Kk, 1, Lambda); % state indices of each chain

% for exchange steps
even=1:2:(W-1);
odd=2:2:(W-1);
pair=[1:(W-1);  2:W]';
evenodd=0;
update_first=zeros(W,1);

% simulate chains in parallel
params = ones(W, 2); % initial (sub-)chain k=1 on each worker
params(:,2) = zeros(W, 1); % initial hold time h=0 on each worker
ti = 1;
%initialise chains
if strcmp(init, 'prior')
    Xin = gamrnd_mixture(alpha, theta, Kk, W);
else
    Xin = init+zeros(Kk, W);
end
 
active_k = zeros(W, 2);

for t = tis
    nold=n;
    X = zeros(Kk, deltat+100, W);
    
    %% PARALLEL LOCAL MOVES %%
    parfor w = 1:W
        par = params;
        update_1 = update_first(w);
        t1 = t; % deadline
        lambda=lambdas(w); % temperature
        X1 = zeros(Kk, deltat+100);
        
        % initialise/update parameters
        ti1 = ti;
        k1 = par(w,1);
        h1 = par(w,2);
        n1 = ones(Kk,1);
        
        % initialise/update state
        x = Xin(:,w);
        X1(:, 1) = x;
        
        %% SIMULATE REAL-TIME MARKOV JUMP PROCESS UNTIL REAL TIME t1 %%
        while ti1 <= t1
            % update k1-th state directly if previous hold time ended
            % at the same time as the deadline
            if(~update_1)
                % HOLD k1-th chain state
                if h1<1 %make sure state isn't already on hold
                    h1 = h1 + gamrnd_mixture(x(k1)^p./theta, theta, 1, 1, 1);
                end
                
                while h1 > 1 && ti1 <= t1
                    ti1 = ti1 + 1;
                    h1 = h1 - 1;
                end
                if(h1 <= 1)
                    update_1=1;
                end
            end
            % UPDATE k1-th chain state
            if (ti1 <= t) %do not update yet if deadline has occured
                n1(k1) = n1(k1) + 1;
                % update k1-th state using random walk Metropolis
                x_star = normrnd(x(k1), rho);
                % odds ratio
                A1 = gampdf_mixture_temp(x_star, lambda, Lambda, alpha, theta);
                A2 = gampdf_mixture_temp(x(k1), lambda, Lambda, alpha, theta);
                A = A1/A2;
                % accept proposal with probability A
                if(unifrnd(0,1) < A)
                    x(k1) = x_star;
                end
                % record update
                X1(k1, n1(k1)) = x(k1);
                update_1=0;
                
                % next chain
                k1 = k1 + 1;
                if (k1 > Kk)
                    k1 = 1;
                end
            end
        end
        active_k(w,:) = [k1 h1];
        update_first(w) = update_1;
        n(:,:,w) = n1;
        X(:,:,w) = X1;
    end
    
    %% PERFORM EXCHANGE MOVES %%
    
    % identify the currently active (sub-)chains k on each worker
    k = active_k(:,1);
    
    % bias correction (or not)
    if(correct) % discard active chains k on each worker
        L = Kk-1;
        exK = zeros(W, L);
        for w=1:W
            exK(w, :) =  [1:(k(w)-1) (k(w)+1):Kk];
        end
    else % keep all chains (no bias correction)
        L = Kk;
        exK = reshape(repelem(1:Kk, W), W, Kk);
    end
    
    % extract current state (and corresponding index) for all chains
    XE = zeros(W, L);
    nE = zeros(W, L);
    
    for w = 1:W
        for l=1:L
            nE(w, l) = n(exK(w,l), :, w);
            XE(w, l) = X(exK(w,l), nE(w, l), w);
        end
    end
    
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
    
    % perform exchange moves across all workers
    for l=1:L  % exchange moves for various (sub-)chains within workers performed in parallel for speed
        XE1 = XE(:,l);
        nE1 = nE(:, l);
        if(~correct) %Note: code not adapted for Kk>1
            hE1 = active_k(:,2);
        end
        lbds=lambdas;
        pair1=pair;
        
        nE1(swaps) = nE1(swaps)+1;
        nE1(swaps+1) = nE1(swaps+1)+1;
        
        % swapping more than one chain at once
        swap = pair1(swaps,:);
        i=swap(:,1);
        j=swap(:,2);
        % odds ratio
        A1=gampdf_mixture_temp(XE1(i)', lbds(j), Lambda, alpha, theta).*gampdf_mixture_temp(XE1(j)', lbds(i), Lambda, alpha, theta);
        A2=gampdf_mixture_temp(XE1(i)', lbds(i), Lambda, alpha, theta).*gampdf_mixture_temp(XE1(j)', lbds(j), Lambda, alpha, theta);
        A=A1./A2;
        % accept or reject swaps
        mkswap = unifrnd(0,1,size(A)) < A;
        % perform the accepted swaps
        x2 = XE1(swap(mkswap,2));
        XE1(swap(mkswap,2)) = XE1(swap(mkswap,1));
        XE1(swap(mkswap,1)) = x2;
        if(~correct) % also swap hold times if no bias correction
            hE1(swap(mkswap,2)) = hE1(swap(mkswap,1));
            hE1(swap(mkswap,1)) = h2;
            active_k(:,2) = hE1;
        end
        
        XE(:,l) = XE1;
        nE(:,l) = nE1;
    end
    
    params=active_k;
    % record swaps
    for w = 1:W
        for l=1:L
            n(exK(w,l), :, w) = nE(w, l);
            X(exK(w,l), nE(w, l), w) = XE(w, l);
        end
        
        for ll=1:Kk
            Xin(ll,w) = X(ll, n(ll, :, w), w);
        end
    end
    
    % record updates and swaps that occured in the time interval
    nnew = n;
    n = n+nold-1;
    for w = 1:W
        for ll=1:Kk
            Xout(ll, nold(ll, :, w):n(ll, :, w), w) = X(ll, 1:nnew(ll, :, w), w);
        end
    end
    
    % prepare to resume update
    ti=t;
    
end
nmax = max(n(:));
Xout = Xout(:,1:nmax,:);
end
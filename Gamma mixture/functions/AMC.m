function [X, Ac, n] = AMC(K, T, N, C, rho, alpha, theta, p, init)
% Standard (Anytime) Parallel Monte Carlo algorithm
% Inputs:   K: # chains
%           T: units of virtual time for which to run algorithm
%           N:  # slots allocated to matrix X to record samples
%           C:  run the algorithm C times in parallel
%           rho: standard deviation for Metropolis update
%           alpha, theta: parameters of Gamma mixture
%           p: computational complexity
%           init: point of initialisation of chains
% Outputs:  X: all chains
%           Ac: accepted local moves for all chains
%           n: sample size obtained for each chain

if ~exist('init', 'var')
    init='prior';
end

X = zeros(K, N, C);
Ac = zeros(K, N, C);
n = zeros(K, C);

parfor c=1:C
    X1 = zeros(K, N);
    Ac1 = zeros(K, N);    
   
    % initialise K chains
    if strcmp(init, 'prior')
        x = gamrnd_mixture(alpha, theta, K, 1);
    else
        x = zeros(K, 1) + init;
    end
    
    % simulate chains one by one
    k = 1;
    t = 1;
    
    h = 0;
    n1 = ones(K, 1);
    %     tic
    %% SIMULATE REAL-TIME MARKOV JUMP PROCESS %%
    while t <= T
        
        % simulate chains one by one
        % hold k-th state
            h = h + gamrnd_mixture(x(k)^p./theta, theta, 1, 1, 1);
        while h > 1 && t <= T
            t = t + 1;
            h = h - 1;
        end
        if t <= T %do not update if deadline has occured
            % update k-th state using random walk Metropolis
            x_star = normrnd(x(k), rho);
            % odds ratio
            A1 = gampdf_mixture(x_star, alpha, theta, 0.5);
            A2 = gampdf_mixture(x(k), alpha, theta, 0.5);
            A = A1/A2;
            % accept proposal with probability A
            if(unifrnd(0,1) < A)
                x(k) = x_star;
                Ac1(k, n1(k))=1;
            end
                     
            X1(k, n1(k)) = x(k);
            
            % next chain
            n1(k) = n1(k) + 1;
            k = k + 1;
            if (k > K)
                k = 1;
            end
        end
        
    end
    X(:,:,c) = X1;
    Ac(:,:,c) = Ac1;
    n(:,c) = n1;
end
n = n-1;
maxn = max(n(:));
X = X(:,1:maxn,:);
Ac = Ac(:,1:maxn,:);
end
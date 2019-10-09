function [Theta, X, Rej, n, ne, nsw, sw, TM] = LV_ABC_standard_exchange_1(params, exchange, LV, observations)
% Single Processor Standard ABC with exchange moves for Lotka-Volterra model
% Inputs:       K: # chains
%               LV: Lotka-Volterra model settings
%               observations: observations and simulation settings
%               params: various parameters (real time schedule, prior,
%               epsilon)
%               exchange: exchange move parameters
% Outputs:      Theta, X: cold chain(s) and corresponding simulations
%               n: sample size obtained for each chain 
%               ne: first sample after burn in for each chain 
%               nsw: total # exchange moves performed   
%               sw: total # exchange moves accepted 
%               TM: timeline of local and exchange moves

nET = observations.ET-1; 
 

tis = exchange.deltat:(exchange.deltat):params.N;
iter = 1:1000:params.N;

% initialise various variables
ti = 1; k = 1; tm = 0; j = 1; J = length(tis); 
nsw = 0; sw = 0; burning = 1; TM = zeros(J,4); tidx = 0;
theta = params.theta_in; % initialise theta for each chain
n = ones(params.K,1); ne = n;
Theta = zeros(params.K, 3, params.N); Theta(:, :, 1) = theta;
X = zeros(params.K, nET, params.N); Rej = zeros(params.K, params.N); 

% more exchange parameters
even = 1:2:(params.K-1); sze = size(even, 2);
odd = 2:2:(params.K-1); szo = size(odd, 2);
evenodd = 1;

deadline.T = params.T; deadline.all = tic;
% initialise x
x = get_x_multi(params.K, LV, observations, theta, deadline); X(:,:, 1) = x;

while (j <= J)&&(tm<params.T)
    t = tis(j); j = j+1; tidx = tidx+1; 
    %% LOCAL MOVES %%  
    a = tic;
    while (ti <= t)&&(tm<params.T) %local move
        if(sum(ti==iter)~=0)
            fprintf('sample %d, time %f \n', n(1), toc(deadline.all))
        end
        n(k) = n(k) + 1;
        [theta(k,:), x(k,:), rej] = LV_local_moves_standard(deadline, theta(k,:), x(k,:), observations, params, LV, params.epsilon(k), params.SIGMA(:,k));
        if(rej==-2)
            % retry current move with another proposal if taking too long
            n(k) = n(k) - 1;
        else
            % update kth chain
            Theta(k, :, n(k)) = theta(k,:);
            X(k, :, n(k)) = x(k,:);
            Rej(k, n(k)) = rej;
            % next time step
            ti = ti + (rej~=2); % do not count rejection due to prior as the next step
        end
                     
        % next chain        
        if(rej~=-2) % unless chain reset because taking too long
            k = k + 1;
            if (k > params.K)
                k = 1;
            end
        end
        tm = toc(deadline.all);
    end
    TM(tidx,:) = [toc(a) 1 tm n(1)];
    
    %%%%%%%%%%%%%%%%%%%%%
    % keeping track of time and burn-in
    if burning&&(tm>params.burnin)
        ne = n;
        burning = 0;
    end
    %%%%%%%%%%%%%%%%%%%%%
    
    %% EXCHANGE MOVES %%
    % perform exchange moves before resuming update
    if(tm<params.T)
        tidx = tidx + 1;
        b = tic;        
        % select chains to attempt to swap
        if evenodd
            swap = exchange.pair(even,:);
            ns = sze;
            evenodd = 0;
        else
            swap = exchange.pair(odd,:);
            ns = szo;
            evenodd = 1;
        end
        
        % perform exchange moves
        nsw = nsw + ns;
        [theta, x, rej_swap, n, swap] = LV_exchange_moves_m(theta, x, n, observations, params.epsilon, ns, swap);
        sw = sw + sum(rej_swap==4);
        
        % update chains
        for i=1:ns
            Rej(swap(i, 1), n(swap(i, 1))) = rej_swap(i);
            Rej(swap(i, 2), n(swap(i, 2))) = rej_swap(i);
            Theta(swap(i, 1), :, n(swap(i, 1))) = theta(swap(i, 1),:);
            Theta(swap(i, 2), :, n(swap(i, 2))) = theta(swap(i, 2),:);
            X(swap(i, 1), :, n(swap(i, 1))) = x(swap(i, 1), :);
            X(swap(i, 2), :, n(swap(i, 2))) = x(swap(i, 2), :);
        end
        TM(tidx,:) = [toc(b) 2 tm n(1)];
    end
end
nmax = max(n(:));
Theta = Theta(:,:,1:nmax); X = X(:,:,1:nmax);
Rej = Rej(:, 1:nmax);TM = TM(TM(:,1)~=0,:);
end
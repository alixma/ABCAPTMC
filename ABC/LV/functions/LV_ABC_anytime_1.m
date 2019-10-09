function [Theta, X, Rej, n, ne, nsw, sw, TM] = LV_ABC_anytime_1(K, observations, LV, params, exchange)
% Single Processor Anytime ABC with exchange moves for Lotka-Volterra model
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

% initialise theta for all chains
theta = params.theta_in; % exprnd(repmat(1./params.prior, K, 1)); %
% initialise various other variables
Theta = zeros(K, 3, params.N); Theta(:, :, 1) = theta;
X = zeros(K, nET, params.N); Rej = zeros(K, params.N); 
TM = zeros(params.N,4); tidx = 0;
n = ones(K, 1); ne = ones(K, 1); nsw = 0; sw = 0;
k = 1; resume = 0; resumesim.resume = 0; t = 0;
rescount = zeros(K, 1); burning = 1; 

% more exchange parameters
even = 1:2:(params.K-2); sze = size(even, 2);
odd = 2:2:(params.K-2); szo = size(odd, 2);
evenodd = 1;

% start the clock
deadline.T = params.T; deadline.deltat = exchange.deltat;
deadline.all = tic;
% initialise x
x = get_x_multi(K, LV, observations, theta, deadline); X(:,:, 1) = x;

while t<=params.T
    ti=0; tidx = tidx + 1; r = tic;       
    deadline.r = r; deadline.ti = ti; deadline.t = t;    
    %% LOCAL MOVES %%
    while(ti<=exchange.deltat)&&(t<=params.T)
        if(~resume)
            params.theta_current = theta(k,:);
            params.x_current = x(k,:);
        end
        % perform local moves for deltat seconds
        [params.theta_current, params.x_current, resume, resumesim, ti, rej] = LV_local_moves_anytime(deadline, observations, params, LV, exchange, params.epsilon(k), params.SIGMA(:,k), resume, resumesim);        
        
        %%%%%%%%%%%%%%%%%%%%%
        rescount(k) = rescount(k) + resume;
        t = toc(deadline.all); % keeping track of time
        %%%%%%%%%%%%%%%%%%%%%%
        
        if(~resume)
            n(k) = n(k)+1;
            theta(k,:) = params.theta_current;
            x(k,:) = params.x_current;
            Theta(k, :, n(k)) = params.theta_current;
            X(k, :, n(k)) = params.x_current;
            Rej(k, n(k)) = rej;
            
            % switch to next chain
            rescount(k) = 0; % reset count
            k = k+1;
            if(k>K)
                k = 1;
            end
        elseif(rescount(k)>0)
            % if same chain resumes too many times in a row, try again with new proposal
            resume = 0; resumesim.resume = 0;
            rescount(k) = 0;   % reset count
            % fprintf("\n resetting chain %d", k)
        end
        %%%%%%%%%%%%%%%%%%%%%%
    end
    TM(tidx,:) = [ti 1 t n(1)];
    
    %%%%%%%%%%%%%%%%%%%%%    
    t = toc(deadline.all); % keeping track of time and burn-in
    if burning&&(t>params.burnin)
        ne = n;
        burning = 0;
    end
    %%%%%%%%%%%%%%%%%%%%%%
    
    %% EXCHANGE MOVES %%
    if(t<params.T)
        tidx = tidx + 1;
        b = tic;
        % discard working chain (bias correction)
        exK = [1:(k-1) (k+1):K];        
        thetaE = theta(exK,:); xE = x(exK,:);        
        nE = n(exK); epsE = params.epsilon(exK);        
        % select chains to attempt to swap
        if evenodd
            swap = exchange.pair(even,:);
            ns = sze;
            evenodd=0;
        else
            swap = exchange.pair(odd,:);
            ns = szo;
            evenodd = 1;
        end
        
        % perform exchange moves
        nsw = nsw + ns;
        [thetaE, xE, rej_swap, nE, swap] = LV_exchange_moves_m(thetaE, xE, nE, observations, epsE, ns, swap);
        sw = sw + sum(rej_swap==4);
        
        % record swap        
        theta(exK,:) = thetaE; x(exK,:) = xE; n(exK) = nE;                
        % update chains
        for i=1:ns
            Theta(exK(swap(i, 1)),:, nE(swap(i, 1))) = thetaE(swap(i, 1),:);
            Theta(exK(swap(i, 2)),:, nE(swap(i, 2))) = thetaE(swap(i, 2),:);
            X(exK(swap(i, 1)), :, nE(swap(i, 1))) = xE(swap(i, 1),:);
            X(exK(swap(i, 2)),:, nE(swap(i, 2))) = xE(swap(i, 2),:);
            Rej(exK(swap(i, 1)), nE(swap(i, 1))) = rej_swap(i);
            Rej(exK(swap(i, 2)), nE(swap(i, 2))) = rej_swap(i);
        end        
        TM(tidx,:) = [toc(b) 2 t n(1)];
    end
    
end
nmax = max(n(:));
Theta = Theta(:,:, 1:nmax); X = X(:, :, 1:nmax);
Rej = Rej(:,1:nmax); TM = TM(TM(:,1)~=0,:);
end
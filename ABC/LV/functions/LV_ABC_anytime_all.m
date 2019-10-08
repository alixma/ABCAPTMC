function [Theta, X, Rej, n, ne, nsw, sw, TM, alln] = LV_ABC_anytime_all(K, observations, LV, params, exchange, nocold)
% Single Processor Anytime LV algorithm
% Note:         Outputs all chains
% Inputs:       K: # chains
%               LV: Lotka-Volterra model settings
%               observations: observations and simulation settings
%               params: various parameters (real time schedule, prior,
%               epsilon)
%               exchange: exchange move parameters
% Outputs:      Theta, X: cold chain(s) and corresponding simulations
%               n: sample size obtained for each chain 
%               ne: first sample at minimum epsilon for each chain 
%               nsw: total # exchange moves performed   
%               sw: total # exchange moves accepted   
if(~exist('nocold', 'var'))
    nocold=0;
end

epsilon = params.epsilon;
nET = observations.ET-1;

%initialise theta for all chains
theta = params.theta_in; %params.S_rej(randsample(1:length(params.S_rej), K),:); %exprnd(repmat(1./params.prior, K, 1));

Theta = zeros(K, 3, params.N);
X = zeros(K, nET, params.N);
Rej = zeros(K, params.N);
Theta(:, :, 1) = theta;
TM = zeros(params.N,4); tidx=0;
alln = zeros(params.N,K);
n = ones(K, 1);
ne = ones(K, 1);
nsw=0; sw=0;
burnin=params.burnin;
burning=1;
k=1; resume=0; resumesim.resume=0; t=0;
% ADDED 19/09/19 keep count of how many times a single chain resumes local
% moves after exchange moves
rescount=zeros(K, 1);
%start the clock
notresuming=0;
deadline.T = params.T;

even = 1:2:(params.K-2); sze = size(even, 2);
odd = 2:2:(params.K-2); szo = size(odd, 2);
evenodd = 1;

deadline.deltat = exchange.deltat;
deadline.all=tic;
x = get_x_multi(K, LV, observations, theta, deadline);
X(:,:, 1) = x;



while t<=params.T
    ti=0;
    tidx = tidx+1;
    r = tic; 
    deadline.r = r; deadline.ti=ti;
    deadline.t = t;
    %% LOCAL MOVES %%
    while(ti<=exchange.deltat)&&(t<=params.T)
        if(~resume)
            params.theta_current = theta(k,:);
            params.x_current = x(k,:);
        end
        %perform local moves for deltat seconds
        if(k~=1)||((k==1)&&(~nocold))            
            [params.theta_current, params.x_current, resume, resumesim, ti, rej] = LV_local_moves_anytime(deadline, observations, params, LV, exchange, epsilon(k), params.SIGMA(:,k), resume, resumesim);
        end
        rescount(k) = rescount(k) + resume;
        
        %%%%%%%%%%%%%%%%%%%%%
        t=toc(deadline.all); %keeping track of time
        %%%%%%%%%%%%%%%%%%%%%%
        
        if(~resume)            
            if(k~=1)||((k==1)&&(~nocold))
                notresuming = notresuming+1;
                n(k) = n(k)+1;
                theta(k,:) = params.theta_current;
                x(k,:) = params.x_current;
                Theta(k, :, n(k)) = params.theta_current;
                X(k, :, n(k)) = params.x_current;
                Rej(k, n(k)) = rej;
            end
            % switch to next chain
            rescount(k) = 0; % reset count
            k = k+1;
            if(k>K)
                k=1;
            end
        elseif(rescount(k)>2)
            % if same chain resumes more than twice in a row, get new proposal
            params.theta_current = params.theta_current(1, :);
            params.x_current = params.x_current(1, :);
            resume = 0; resumesim.resume=0;
            rescount(k) = 0;   % reset count
            %fprintf("\n resetting chain %d", k)
        end
        %%%%%%%%%%%%%%%%%%%%%%
    end
    TM(tidx,:) = [ti 1 t n(1)];
    alln(tidx,:) = n;
    %%%%%%%%%%%%%%%%%%%%%    
    t=toc(deadline.all); %keeping track of time
    if burning&&(t>burnin)
        ne = n;
        burning=0;
    end
    %%%%%%%%%%%%%%%%%%%%%%
    %theta
    %% EXCHANGE MOVES %%
    if(t<params.T)
        tidx = tidx+1;
        b=tic;
        % discard working chain (bias correction)
        exK = [1:(k-1) (k+1):K];        
        thetaE = theta(exK,:);
        xE = x(exK,:);
        nE = n(exK);
        epsE = epsilon(exK);
        if evenodd
            swap = exchange.pair(even,:);
            ns = sze;
            evenodd=0;
        else
            swap = exchange.pair(odd,:);
            ns = szo;
            evenodd = 1;
        end
        nsw = nsw+ns;
        [thetaE, xE, rej_swap, nE, swap]=LV_exchange_moves_standard_m(thetaE, xE, nE, observations, epsE, ns, swap);
        sw = sw + sum(rej_swap==4);
        
        %record swap        
        theta(exK,:) = thetaE;
        x(exK,:) = xE;
        n(exK) = nE;
        %update chains
        for i=1:ns
            Theta(exK(swap(i, 1)),:, nE(swap(i, 1))) = thetaE(swap(i, 1),:);
            Theta(exK(swap(i, 2)),:, nE(swap(i, 2))) = thetaE(swap(i, 2),:);
            X(exK(swap(i, 1)), :, nE(swap(i, 1))) = xE(swap(i, 1),:);
            X(exK(swap(i, 2)),:, nE(swap(i, 2))) = xE(swap(i, 2),:);
            Rej(exK(swap(i, 1)), nE(swap(i, 1))) = rej_swap(i);
            Rej(exK(swap(i, 2)), nE(swap(i, 2))) = rej_swap(i);
        end
        
        TM(tidx,:) = [toc(b) 2 t n(1)];
        alln(tidx,:) = n;
    end
    
end
nmax = max(n(:));
Theta = Theta(:,:, 1:nmax);
X = X(:, :, 1:nmax);
Rej = Rej(:,1:nmax);
TM = TM(TM(:,1)~=0,:);
end
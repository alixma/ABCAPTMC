function [Theta, X, Rej, n, ne, nsw, sw, TM, epsilon] = LV_ABC_anytime_all_adaptive(K, observations, LV, params, exchange, nocold)
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
TM = zeros(params.N,2); tidx=0;
alln = zeros(params.N,K);
n = ones(K, 1);
ne = ones(K, 1); n_idx = n; d=0.034;
nsw=0; sw=0;
burnin=params.burnin;
burning=1;
k=1; resume=0; resumesim.resume=0; t=0;
%start the clock
notresuming=0;
deadline.T = params.T;

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
            k = k+1;
            if(k>K)
                k=1;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%
%         if(~record) % record first n with minimum epsilon
%             if(sum(x==1e3))
%                 ne = n;
%                 record=1;
%             end
%         end
        %%%%%%%%%%%%%%%%%%%%%%
    end
    TM(tidx,:) = [ti 1];
    alln(tidx,:) = n;
    %%%%%%%%%%%%%%%%%%%%%    
    t=toc(deadline.all); %keeping track of time
    if burning&&(t>burnin)
        ne = n;
        burning=0;
    end
    
    % check acceptance rate of local moves
    if sum(n(2:end)-n_idx(2:end)>=50)==K-1 %if all warm chains have 50 more samples than last time
        for cidx=2:K
            Rej_idx = Rej(cidx,1:n(cidx));
            a_idx = sum(Rej_idx==0)/sum(Rej_idx<3);
            fprintf('\n Chain %d: acceptance rate = %f --> ', cidx, a_idx)
            % update if not between 20%-26%
            if (a_idx<0.234-d)
                % too low, increase epsilon                
                epsilon(cidx) = epsilon(cidx) + abs(normrnd(0, params.SIGMA(1,cidx)));
                fprintf('increasing epsilon to %f \n', epsilon(cidx))
            elseif (a_idx>0.234+d)
                % too high, reduce epsilon
                epsilon(cidx) = epsilon(cidx) - abs(normrnd(0, params.SIGMA(1,cidx)));
                fprintf('reducing epsilon to %f \n', epsilon(cidx))
            else
                fprintf('all good \n')
            end
        end
        n_idx = n;
        params.epsilon = epsilon;
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
        nsw = nsw+1;
        % randomly select pair to attempt swap
        swap = sort(exchange.pair(randsample(exchange.ipairs, 1),:));
        [thetaE, xE, rej_swap, nE, swap]=LV_exchange_moves_clock(deadline, thetaE, xE, nE, observations, LV, exchange, epsE, swap);
        sw = sw + 1*(rej_swap==4);
        
        %record swap        
        theta(exK,:) = thetaE;
        x(exK,:) = xE;
        n(exK) = nE;
        %update chains
        Theta(exK(swap(1)),:, nE(swap(1))) = thetaE(swap(1),:);
        Theta(exK(swap(2)),:, nE(swap(2))) = thetaE(swap(2),:);
        X(exK(swap(1)), :, nE(swap(1))) = xE(swap(1),:);
        X(exK(swap(2)),:, nE(swap(2))) = xE(swap(2),:);
        Rej(exK(swap(1)), nE(swap(1))) = rej_swap;
        Rej(exK(swap(2)), nE(swap(2))) = rej_swap;
        
        TM(tidx,:) = [toc(b) 2];
        alln(tidx,:) = n;
    end
    
end
nmax = max(n(:));
Theta = Theta(:,:, 1:nmax);
X = X(:, :, 1:nmax);
Rej = Rej(:,1:nmax);
TM = TM(TM(:,1)~=0,:);
end
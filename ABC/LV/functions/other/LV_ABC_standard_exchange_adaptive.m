function [Theta, X, ne, n, Rej, nsw, sw, TM, epsilon] = LV_ABC_standard_exchange_adaptive(params, exchange, LV, observations, nocold)
% Standard ABC algorithm for Lotka-Volterra example

if(~exist('nocold', 'var'))
    nocold=0;
end

nET = observations.ET-1; y = observations.y; epsilon = params.epsilon;
T = params.T; deadline.T = params.T; deltat = exchange.deltat; N=params.N; K = params.K; %prior = params.prior;
 burning=1;

tis = (1+deltat):deltat:N;
iter = 1:1000:N;

%exchange parameters
ipairs = exchange.ipairs; pair = exchange.pair;

% initialise theta and x for each chain
ti=1; k=1; tm=0;
j=1; J = length(tis); 
nsw = 0; sw = 0;
theta = params.theta_in; %params.S_rej(randsample(1:length(params.S_rej), K),:);
TM = zeros(J,2); tidx=0;
n = ones(K,1); ne=n; n_idx = n; d=0.034;
Theta = zeros(K, 3, N);
Rej = zeros(K, N);
X = zeros(K, nET, N);
Theta(:, :, 1) = theta;


deadline.all=tic;
x = get_x_multi(K, LV, observations, theta, deadline);
X(:,:, 1) = x;

while (j <= J)&&(tm<T)
    t = tis(j);
    j = j+1;
    tidx = tidx+1;
    %tm=toc;
    %     if(mod(t-1, T{1}/20)==0)
    %         fprintf('iteration %d of %d \n', t, T{1})
    %     end
    %% LOCAL MOVES %%  
    a=tic;
    while (ti <= t)&&(tm<T) %local move
        if(sum(ti==iter)~=0)
            fprintf('sample %d, time %f \n', n(1), toc(deadline.all))
        end
        if(k~=1)||((k==1)&&(~nocold))
            n(k) = n(k) + 1; 
            [theta(k,:), x(k,:), rej]=LV_local_moves_standard(deadline, theta(k,:), x(k,:), observations, params, LV, epsilon(k), params.SIGMA(:,k));
            Theta(k, :, n(k)) = theta(k,:);
            X(k, :, n(k)) = x(k,:);
            Rej(k, n(k)) = rej;
            % next time step
            ti = ti+(rej~=2);%1; %do not count rejection due to prior a the next step
        end
             
        % next chain
        k = k+1;
        if (k > K)
            k = 1;
        end
        tm=toc(deadline.all);
    end
    TM(tidx,:) = [toc(a) 1];
    %t=toc(deadline.all); %keeping track of time
    if burning&&(t>params.burnin)
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
    
    
    
    %% EXCHANGE MOVES %%
    % perform exchange moves before resuming update
    if(tm<T)
        tidx = tidx+1;
        b=tic;
        nsw=nsw+1;
        [theta, x, rej_swap, n, swap]=LV_exchange_moves_standard(deadline, theta, x, n, observations, LV, exchange, epsilon);
        sw = sw + 1*(rej_swap==4);
        
        % record swap
        Rej(swap(1), n(swap(1))) = rej_swap;
        Rej(swap(2), n(swap(2))) = rej_swap;
        Theta(swap(1), :, n(swap(1))) = theta(swap(1),:);
        Theta(swap(2), :, n(swap(2))) = theta(swap(2),:);
        X(swap(1), :, n(swap(1))) = x(swap(1), :);
        X(swap(2), :, n(swap(2))) = x(swap(2), :);
        TM(tidx,:) = [toc(b) 2];
    end
end
nmax = max(n(:));
Theta = Theta(:,:,1:nmax);
X = X(:,:,1:nmax);
Rej = Rej(:, 1:nmax);
TM = TM(TM(:,1)~=0,:);
end
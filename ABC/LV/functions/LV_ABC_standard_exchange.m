function [Theta, X, ne, n, Rej, nsw, sw, TM] = LV_ABC_standard_exchange(params, exchange, LV, observations, nocold)
% Standard ABC algorithm for Lotka-Volterra example

if(~exist('nocold', 'var'))
    nocold=0;
end

nET = observations.ET-1; y = observations.y; epsilon = params.epsilon;
T = params.T; deadline.T = params.T; deltat = exchange.deltat; N=params.N; K = params.K; %prior = params.prior;
 burning=1;

tis = (deltat):(deltat):N;
iter = 1:1000:N;

%exchange parameters
ipairs = exchange.ipairs; pair = exchange.pair;

% initialise theta and x for each chain
ti=1; k=1; tm=0;
j=1; J = length(tis); 
nsw = 0; sw = 0;
theta = params.theta_in; %params.S_rej(randsample(1:length(params.S_rej), K),:);
TM = zeros(J,4); tidx=0;
n = ones(K,1); ne=n;
Theta = zeros(K, 3, N);
Rej = zeros(K, N);
X = zeros(K, nET, N);
Theta(:, :, 1) = theta;

%swap multiple pairs at once
even=1:2:(K-1); sze = size(even, 2);
odd=2:2:(K-1); szo = size(odd, 2);
evenodd=1;

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
            if(rej==-2)
                % reset
                n(k) = n(k) - 1; 
            else
                % update kth chain
                Theta(k, :, n(k)) = theta(k,:);
                X(k, :, n(k)) = x(k,:);
                Rej(k, n(k)) = rej;
                % next time step
                ti = ti + (rej~=2); %do not count rejection due to prior or resetting as the next step 
            end
        end
             
        % next chain
        % unless chain reset because taking too long
        if(rej~=-2)
            k = k+1;
            if (k > K)
                k = 1;
            end
        end
        tm=toc(deadline.all);
    end
    TM(tidx,:) = [toc(a) 1 tm n(1)];
    %t=toc(deadline.all); %keeping track of time
    if burning&&(tm>params.burnin)
        ne = n;
        burning=0;
    end
    
    %% EXCHANGE MOVES %%
    % perform exchange moves before resuming update
    if(tm<T)
        tidx = tidx+1;
        b=tic;        
        if evenodd
            swap = exchange.pair(even,:);
            ns = sze;
            evenodd=0;
        else
            swap = exchange.pair(odd,:);
            ns = szo;
            evenodd = 1;
        end
        nsw=nsw+ns;
        [theta, x, rej_swap, n, swap]=LV_exchange_moves_standard_m(theta, x, n, observations, epsilon, ns, swap);
        sw = sw + sum(rej_swap==4);
        
        % record swaps
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
Theta = Theta(:,:,1:nmax);
X = X(:,:,1:nmax);
Rej = Rej(:, 1:nmax);
TM = TM(TM(:,1)~=0,:);
end
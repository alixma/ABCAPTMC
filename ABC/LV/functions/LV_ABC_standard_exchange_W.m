function [Theta_out, X_out, n_out, ne_out, Rej_out, nsw, sw, TM] = LV_ABC_standard_exchange_W(params, exchange, LV, observations)
% Multiple Processor Standard ABC algorithm for Lotka-Volterra example
% Inputs:       params: various parameters (real time schedule, prior,
%               epsilon)
%               exchange: exchange move parameters
%               LV: Lotka-Volterra model settings
%               observations: observations and simulation settings
% Outputs:      Theta, X: cold chain(s) and corresponding simulations
%               n: sample size obtained for each chain 
%               ne: first sample at minimum epsilon for each chain 
%               nsw: total # between worker exchange moves performed   
%               sw: total # between worker exchange moves accepted  
%               TM: times spent performing local moves (overall and on each provessor) and exchange moves

nET = observations.ET-1; 

reverseStr = '';

% change epsilon and sigma to an array
epsilon = reshape(params.epsilon, exchange.Kk, params.W);
SIGMA = reshape(params.SIGMA, 3, exchange.Kk, params.W);

tis = (1+exchange.deltat):exchange.deltat:params.N;

% exchange parameters
% for swaps within workers
pair_w = [1:(exchange.Kk-1); 2:(exchange.Kk)]';
even_w = 1:2:(exchange.Kk-1); sze_w = size(even_w, 2);
odd_w = 2:2:(exchange.Kk-1); szo_w = size(odd_w, 2);

% for swaps between workers
pair_b = [1:(params.K-1);  2:params.K]';
even_b = 1:2:(params.K-1); sze_b = size(even_b, 2);
twr_even = repmat(1:ceil(params.K/2), 2, 1); twr_even=twr_even(:);
odd_b = 2:2:(params.K-1); szo_b = size(odd_b, 2);
twr_odd = repmat(1:ceil(params.K/2), 2, 1); twr_odd=[0; twr_odd(:)];
tww = params.tww; twk = params.twk; evenodd=1;

% initialise theta and x for each chain
ti = ones(1, params.W); nsw = 0; sw = 0;
TM = zeros(params.N, params.W+2); tidx=1;
theta = zeros(exchange.Kk, 3, params.W);
x = zeros(exchange.Kk, nET, params.W);
deadline_init.T = 100; deadline.T = params.T;
deadline_init.all = tic;
for kk=1:exchange.Kk    
    th = params.S_rej(randsample(1:length(params.S_rej), params.W),:); %exprnd(repmat(1./prior, params.W, 1));
    theta(kk,:,:) = th';
    x(kk,:,:) = get_x_multi(params.W, LV, observations, th, deadline_init)';
end
clear('th', 'deadline_init');

% tidier for output
Theta_out = zeros(params.K, 3, params.N);
Rej_out = zeros(params.K, params.N);
X_out = zeros(params.K, nET, params.N);

% less tidy for computations
n = zeros(exchange.Kk, params.W); ne = n;
Theta = zeros(exchange.Kk, 3, params.N, params.W); Theta(:, :, 1, :) = theta;
X = zeros(exchange.Kk, nET, params.N, params.W); X(:, :, 1, :) = x;
Rej = zeros(exchange.Kk, params.N, params.W);
within = ones(1, params.W); 

burning = 1; nsw_w = zeros(1, params.W); sw_w = zeros(1, params.W);
j = ones(1, params.W); tm = 0; J = length(tis); kk = ones(1, params.W);
deadline.all = tic; identifier = [1 2 1]';
exchange.pair_w = pair_w;
    
while (j(1) <= J)&&(tm<deadline.T)
    t_c = zeros(1+2*within(1), params.W);
    
    %% PARALLEL MOVES %%
    exchange.evenodd = evenodd; a = tic;    
    parfor w=1:params.W 
        j1 = j(w); tis1 = tis; t_c1 = t_c(:, w); tidx1 = 1; 
        kk1 = kk(w); exch1 = exchange; sw1 = 0; nsw1 = 1;
        ti1=ti(w); params1 = params; tm1 = tm; deadline1=deadline;
        Theta1 = Theta(:,:,:, w); theta1 = theta(:,:, w);
        X1 = X(:,:,:, w); x1 = x(:,:, w);
        Rej1 = Rej(:,:, w); n1 = n(:,w);
        epsilon1 = epsilon(:, w);
        SIGMA1 = SIGMA(:,:,w);
        within1 = within(w);  something=1;
        while(something)
            something=within1*(tm1<deadline1.T);
            t1 = tis1(j1); j1 = j1+1;
            % fprintf('\n Worker %d: ', w)
            %% LOCAL MOVES %%
            a1 = tic;
            while (ti1 <= t1)&&(tm1<deadline1.T)                 
                n1(kk1) = n1(kk1) + 1;
                [theta1(kk1, :), x1(kk1, :), rej] = LV_local_moves_standard(deadline1, theta1(kk1, :), x1(kk1, :), observations, params1, LV, epsilon1(kk1), SIGMA1(:,kk1));
                Theta1(kk1, :, n1(kk1)) = theta1(kk1, :);
                X1(kk1, :, n1(kk1)) = x1(kk1, :);
                Rej1(kk1, n1(kk1)) = rej;
                
                % next time step
                ti1 = ti1 + (rej~=2); %do not count rejection due to prior a the next step

                % next chain on worker
                kk1 = kk1+1;
                if (kk1 > exchange.Kk)
                    kk1 = 1;
                end
                tm1 = toc(deadline1.all);
            end
            t_c1(tidx1) = toc(a1); tidx1 = tidx1+1;
            
            %% WITHIN WORKER EXCHANGE MOVES %%
            if(within1)&&(tm1<deadline1.T)
                within1=0;
                b1=tic;
                
                if exch1.evenodd
                    swap = exch1.pair_w(even_w,:);
                    sz = sze_w;
                else
                    swap = exch1.pair_w(odd_w,:);
                    sz = szo_w;
                end
                
                % perform exchange moves
                nsw1 = nsw1 + sz
                [theta1, x1, rej_swap, n1, swap_w] = LV_exchange_moves_m(theta1, x1, n1, observations, epsilon1', sz, swap);
                % fprintf('| swap chains %d and %d ', swap_w(:,1), swap_w(:,2));
                sw1 = sw1+sum(rej_swap==4);
                
                % record swap
                for i=1:sz
                    Rej1(swap_w(i, 1), n1(swap_w(i, 1))) = rej_swap(i);
                    Rej1(swap_w(i, 2), n1(swap_w(i, 2))) = rej_swap(i);
                    Theta1(swap_w(i, 1), :, n1(swap_w(i, 1))) = theta1(swap_w(i, 1),:);
                    Theta1(swap_w(i, 2), :, n1(swap_w(i, 2))) = theta1(swap_w(i, 2),:);
                    X1(swap_w(i, 1), :, n1(swap_w(i, 1))) = x1(swap_w(i, 1), :);
                    X1(swap_w(i, 2), :, n1(swap_w(i, 2))) = x1(swap_w(i, 2), :);
                end
                t_c1(tidx1) = toc(b1); tidx1 = tidx1+1;
            end
        end        

        Theta(:,:,:,w) = Theta1; theta(:, :, w) = theta1;
        X(:, :, :, w) = X1; x(:, :, w) = x1;
        Rej(:,:, w) = Rej1; n(:, w) = n1;        
        ti(w) = ti1; kk(w) = kk1;        
        t_c(:, w) = t_c1; nsw_w(w) = nsw1; sw_w(w) = sw1;         
        within(w) = within1; j(w) = j1;        
    end
    t_a = toc(a); sz = size(t_c, 1); tidx_new = tidx + sz-1;
    TM(tidx:tidx_new,:) = [repmat(t_a, sz, 1) t_c identifier(1:sz)];
    tidx = tidx_new + 1; within = ones(1, params.W);
    nsw = nsw + sum(nsw_w); sw = sw + sum(sw_w);
    
    %%%%%%%%%%%%%%%%%%%%%%
    % keeping track of time and burn-in
    tm = toc(deadline.all);
    if burning&&(tm>params.burnin)
        ne = n;
        burning=0;
    end
    TimeRemaining = params.T-tm;
    msg = sprintf('Time remaining: %3.1f seconds', max(0, TimeRemaining));
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    %%%%%%%%%%%%%%%%%%%%%%
    
    %% BETWEEN WORKER EXCHANGE MOVES %%        
    if(tm<deadline.T)        
        % perform exchange moves between workers before resuming update
        % select pairs to attempt swap 
        if evenodd
            evenodd = 0;
            swap_b = pair_b(even_b,:);
            sz = sze_b; twr = twr_even;            
        else
            evenodd = 1;
            swap_b = pair_b(odd_b,:);
            sz = szo_b; twr = twr_odd;            
        end
        nsw = nsw + sz;
        
        % flatten arrays for exchange moves
        theta_swap = [reshape(theta(:, 1, :), params.K, 1), reshape(theta(:, 2, :), params.K, 1), reshape(theta(:, 3, :), params.K, 1)];
        x_swap = zeros(params.K, nET);
        for i=1:nET
            x_swap(:,i) = reshape(x(:, i, :), params.K, 1);
        end
        epsilon_swap = epsilon(:)'; n_swap = n(:);
        
        
        % perform exchange moves
        b = tic;
        [theta_swap, x_swap, rej_swap, n_swap, ~] = LV_exchange_moves_m(theta_swap, x_swap, n_swap, observations, epsilon_swap, sz, swap_b);
        TM(tidx,:) = [toc(b) zeros(1, params.W) 3]; tidx = tidx+1;        
        sw = sw + sum(rej_swap==4);
        
        % redistribute & record
        n = reshape(n_swap, exchange.Kk, params.W);
        for i=swap_b(:)'
            % redistribute
            theta(twk(i), :, tww(i)) = theta_swap(i, :);
            x(twk(i), :, tww(i)) = x_swap(i,:);

            % record
            Rej(twk(i), n_swap(i), tww(i)) = rej_swap(twr(i));
            Theta(twk(i), :, n_swap(i), tww(i)) = theta_swap(i,:);
            X(twk(i), :, n_swap(i), tww(i)) = x_swap(i,:);
        end
                
    end
end
    fprintf('\n')
% tidy up
k = 1;
for w=1:params.W
    for kk1=1:exchange.Kk
        X_out(k,:, :) = X(kk1, :, :, w);
        Theta_out(k, :, :) = Theta(kk1, :, :, w);
        Rej_out(k, :) = Rej(kk1, :, w);
        k = k + 1;
    end
end
ne_out = ne(:);
n_out = n(:); nmax = max(n_out);
Theta_out = Theta_out(:,:,1:nmax);
X_out = X_out(:,:,1:nmax);
Rej_out = Rej_out(:, 1:nmax);
TM = TM(TM(:,params.W+2)~=0, :);
end
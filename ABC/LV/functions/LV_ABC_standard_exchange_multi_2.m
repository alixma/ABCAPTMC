function [Theta_out, X_out, n_out, ne_out, Rej_out, nsw, sw, TM] = LV_ABC_standard_exchange_multi_2(params, exchange, LV, observations, twice)
% Multiple Processor Standard ABC algorithm for Lotka-Volterra example
% ALTERNATIVE EXCHANGE MOVES
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
%               TM: times spent performing local moves (overall and on each provessor) and exchange moves

if(~exist('twice', 'var'))
    twice=0;
    % requires size of epsilon and sigma to be Kk*W/2
end

nET = observations.ET-1; 
N=params.N; W = params.W;
Kk = exchange.Kk; K = params.K; %prior = params.prior;
deadline.T = params.T;
reverseStr = '';


epsilon = params.epsilon;
SIGMA = params.SIGMA;


%change epsilon to an array
epsilon = reshape(epsilon, Kk, W);
SIGMA = reshape(SIGMA, 3, Kk, W);

tww = params.tww; twk = params.twk; 
tis = (1+exchange.deltat):exchange.deltat:N;
iter = 1:1000:N;

%exchange parameters
if(Kk>1)
    % for swaps within workers
    pair_w = [1:(Kk-1); 2:(Kk)]';
    ipairs_w = 1:size(pair_w,1);
    even_w = 1:2:(Kk-1); sze_w = size(even_w, 2);
    odd_w = 2:2:(Kk-1); szo_w = size(odd_w, 2);

    % for swaps between workers
    pair_b = [1:(K-1);  2:K]';   
    ipairs_b = 1:size(pair_b,1);
    even_b = 1:2:(K-1); sze_b = size(even_b, 2);
    twr_even = repmat(1:ceil(K/2), 2, 1); twr_even=twr_even(:);
    odd_b = 2:2:(K-1); szo_b = size(odd_b, 2);
    twr_odd = repmat(1:ceil(K/2), 2, 1); twr_odd=[0; twr_odd(:)];
else
    ipairs_b = 1:(W-1);
    pair_b = [1:(W-1);  2:W]';
end

% initialise theta and x for each chain
ti=ones(1, W);
nsw = 0; sw = 0;
TM = zeros(params.N, W+2); tidx=1;
% TM_i = zeros(params.N, 2, W); tidx_i = zeros(1, W);
theta = zeros(Kk, 3, W);
x = zeros(Kk, nET, W);
deadline_init.T = 100;
deadline_init.all=tic;
for kk=1:Kk    
    th = params.S_rej(randsample(1:length(params.S_rej), W),:); %repmat(params.theta_in, W, 1); %exprnd(repmat(1./prior, W, 1));
    theta(kk,:,:) = th';
    x(kk,:,:) = get_x_multi(W, LV, observations, th, deadline_init)';
end
clear('th', 'deadline_init');

% tidier for output
Theta_out = zeros(K, 3, N);
Rej_out = zeros(K, N);
X_out = zeros(K, nET, N);

% for exchange moves (between workers)
% collect adjacent temperatures from each worker
% Note: two temperatures must be collected for workers that are not first
% nor last
% theta_swap = cell(W, 1); 
% x_swap = cell(W, 1); 
% n_swap = cell(W, 1); 
% epsilon_swap= cell(W, 1);

% less tidy for computations
n = zeros(Kk, W); ne=n;
Theta = zeros(Kk, 3, N, W);
X = zeros(Kk, nET, N, W);
Rej = zeros(Kk, N, W);
X(:, :, 1, :) = x;
Theta(:, :, 1, :) = theta;
within = ones(1, W); evenodd=1;

burning = 1; sw_w = zeros(1, W);
j=ones(1, W); tm=0; J = length(tis); kk=ones(1, W);
deadline.all=tic; identifier = [1 2 1]';
exchange.pair_w = pair_w; exchange.ipairs_w = ipairs_w;
exchange.even_w=even_w; exchange.odd_w=odd_w;
    
while (j(1) <= J)&&(tm<deadline.T)
    t_c = zeros(1+2*within(1), W);
    
    %     t = tis(j);
%     j = j+1;
    %% PARALLEL MOVES %%

    exchange.evenodd=evenodd; 
    a=tic;
    parfor w=1:W 
        j1 = j(w); tis1 = tis; t_c1 = t_c(:, w); tidx1=1; sw1=0;
%         tidx_i1 = tidx_i(w); TM_i1 = TM_i(:, :, w);
        kk1 = kk(w); exch1 = exchange;
        ti1=ti(w); params1 = params; tm1=tm; deadline1=deadline;
        Theta1 = Theta(:,:,:, w); theta1 = theta(:,:, w);
        X1 = X(:,:,:, w); x1 = x(:,:, w);
        Rej1 = Rej(:,:, w); n1 = n(:,w);
        epsilon1 = epsilon(:, w);
        SIGMA1 = SIGMA(:,:,w); 
        within1 = within(w);  something=1;       
        while(something)
            something=within1*(tm1<deadline1.T);
            t1 = tis1(j1); j1 = j1+1; 
%              fprintf('\n Worker %d: ', w)
            %% LOCAL MOVES %%
            a1=tic;
            while (ti1 <= t1)&&(tm1<deadline1.T) %local move
%                 switch kk1
%                     case 1
%                         fprintf('^')
%                     case 2
%                         fprintf('`')
%                     case 3
%                         fprintf('=')
%                     case 4
%                         fprintf('-')
%                     case 5
%                         fprintf('_')
%                 end
                %             cta1 = toc(a1);
                %             if(sum(ti1==iter)~=0)
                %                                 fprintf('Worker %d: sample %d, time = %.3f \n', w, n1(kk1)+1, toc(deadline1.all))
                %             end
                
%                 if(kk1>1)||((kk1==1)&&(mod(w, 2)==0))
                    n1(kk1) = n1(kk1) + 1;
                    [theta1(kk1, :), x1(kk1, :), rej]=LV_local_moves_standard(deadline1, theta1(kk1, :), x1(kk1, :), observations, params1, LV, epsilon1(kk1), SIGMA1(:,kk1));
                    Theta1(kk1, :, n1(kk1)) = theta1(kk1, :);
                    X1(kk1, :, n1(kk1)) = x1(kk1, :);
                    Rej1(kk1, n1(kk1)) = rej;
                    
                    % next time step
                    ti1 = ti1+(rej~=2); %1; %do not count rejection due to prior a the next step
                    % measure time takin for individual temperatures
                    %                 tidx_i1 = tidx_i1+1;
                    %                 TM_i1(tidx_i1, :) = [toc(a1) kk1];
%                 end
                % next chain on worker
                kk1 = kk1+1;
                if (kk1 > Kk)
                    kk1 = 1;
                end
                tm1=toc(deadline1.all);
            end
            t_c1(tidx1) = toc(a1); tidx1 = tidx1+1;
            
            %% WITHIN WORKER EXCHANGE MOVES %%
            if(within1)&&(tm1<deadline1.T)
                within1=0;
                b1=tic;
                
                if exch1.evenodd
                    swap = exch1.pair_w(exch1.even_w,:);
                    sz = sze_w;
                else
                    swap = exch1.pair_w(exch1.odd_w,:);
                    sz = szo_w;
                end
                
                [theta1, x1, rej_swap, n1, swap_w]=LV_exchange_moves_standard_m(theta1, x1, n1, observations, epsilon1', sz, swap);
%                  fprintf('| swap chains %d and %d ', swap_w(:,1), swap_w(:,2)); 
                sw1 = sum(rej_swap==4);
                
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
        
        % for between chain exchange moves, 
        % (only need the first or last idle chain)
%         switch w
%             case W %if worker index is last pick first
%                 theta_swap{w} = theta1(1,:);
%                 x_swap{w} = x1(1,:);
%                 n_swap{w} =  n1(1);
%                 epsilon_swap{w} = epsilon1(1);
%             case 1 %if worker index is first pick last
%                 theta_swap{w} = theta1(end,:);
%                 x_swap{w} = x1(end,:);
%                 n_swap{w} =  n1(end);
%                 epsilon_swap{w} = epsilon1(end);
%             otherwise %otherwise pick first and last
%                 theta_swap{w} = [theta1(1,:); theta1(end,:)];
%                 x_swap{w} = [x1(1,:); x1(end,:)];
%                 n_swap{w} =  [n1(1); n1(end)];
%                 epsilon_swap{w} = [epsilon1(1); epsilon1(end)];
%         end
%         
        Theta(:,:,:,w) = Theta1; theta(:, :, w) = theta1;
        X(:, :, :, w) = X1; x(:, :, w) = x1;
        Rej(:,:, w) = Rej1;
        n(:, w) = n1;
        ti(w) = ti1;
        kk(w) = kk1;
        t_c(:, w) = t_c1;
        sw_w(w) = sw1; 
        %         TM_i(:, :, w) = TM_i1; tidx_i(w) = tidx_i1;
        within(w) = within1;
        j(w) = j1;
    end
    t_a = toc(a); sz = size(t_c, 1); tidx_new = tidx+sz-1;
    TM(tidx:tidx_new,:) = [repmat(t_a, sz, 1) t_c identifier(1:sz)];
    tidx = tidx_new+1; within = ones(1, W);
    nsw = nsw+W; sw = sw+sum(sw_w);
    
    tm=toc(deadline.all);
    if burning&&(tm>params.burnin)
        ne = n;
        burning=0;
    end
    %% Info dump
    TimeRemaining = params.T-tm;
    msg = sprintf('Time remaining: %3.1f seconds', max(0, TimeRemaining));
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    %%
    
    %% BETWEEN WORKER EXCHANGE MOVES %%        
    if(tm<deadline.T)
        %exchange.pair = pair_b; exchange.ipairs = ipairs_b;
        
        % perform exchange moves between workers before resuming update        
        
        %swap_b = sort(pair_b(randsample(ipairs_b, 1), :));
        if evenodd
            evenodd=0;
            swap_b = pair_b(even_b,:);
            sz = sze_b;
            twr = twr_even;
        else
            evenodd=1;
            swap_b = pair_b(odd_b,:);
            sz = szo_b;
            twr = twr_odd;
        end
        nsw = nsw+sz;
        
        %combine params in one
        theta_swap = [reshape(theta(:, 1, :), K, 1), reshape(theta(:, 2, :), K, 1), reshape(theta(:, 3, :), K, 1)];
        x_swap = zeros(K, nET);
        for i=1:nET
            x_swap(:,i) = reshape(x(:, i, :), K, 1);
        end
        epsilon_swap = epsilon(:)';
        n_swap = n(:);
        % EXCHANGE MOVES
        b=tic;
        [theta_swap, x_swap, rej_swap, n_swap, ~]=LV_exchange_moves_standard_m(theta_swap, x_swap, n_swap, observations, epsilon_swap, sz, swap_b);
        %                      fprintf('\n Exchange moves between workers %d and %d \n', swap_b(i, 1), swap_b(i, 2));
        %                      fprintf('\n Epsilons = %.3f and %.3f \n', epsilon_swap_b(1), epsilon_swap_b(2));
        TM(tidx,:) = [toc(b) zeros(1, W) 3];
        tidx = tidx+1;
        sw = sw+sum(rej_swap==4);
        
        %redistribute & record
        n = reshape(n_swap, Kk, W);
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
%tidy up
k=1;
for w=1:W
    for kk1=1:Kk
        X_out(k,:, :) = X(kk1, :, :, w);
        Theta_out(k, :, :) = Theta(kk1, :, :, w);
        Rej_out(k, :) = Rej(kk1, :, w);
        k=k+1;
    end
end
ne_out = ne(:);
n_out = n(:); nmax = max(n_out);
Theta_out = Theta_out(:,:,1:nmax);
X_out = X_out(:,:,1:nmax);
Rej_out = Rej_out(:, 1:nmax);
TM = TM(TM(:,W+2)~=0, :);
% TM_i = TM_i(1:max(sum(TM_i(:,2, :)~=0)), :, :);
end
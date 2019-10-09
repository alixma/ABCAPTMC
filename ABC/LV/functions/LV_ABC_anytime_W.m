function [Theta_out, X_out, Rej_out, n_out, ne_out, nsw, sw, TM] = LV_ABC_anytime_W(params, exchange, LV, observations)
% Multiple Processor Anytime LV algorithm
% Inputs:       params: various parameters (real time schedule, prior,
%               epsilon)
%               exchange: exchange move parameters
%               LV: Lotka-Volterra model settings
%               observations: observations and simulation settings
% Outputs:      Theta, X: cold chain(s) and corresponding simulations
%               n: sample size obtained for each chain 
%               ne: first sample at minimum epsilon for each chain 
%               nsw: total # exchange moves performed   
%               sw: total # exchange moves accepted  
%               TM: times spent performing local moves (overall and on each provessor) and exchange moves

nET = observations.ET-1; 
reverseStr = '';

% change epsilon and SIGMA to an array
epsilon = reshape(params.epsilon, exchange.Kk, params.W);
SIGMA = reshape(params.SIGMA, 3, exchange.Kk, params.W);

% exchange parameters
% for swaps within workers
pair_w = [1:(exchange.Kk-1); 2:(exchange.Kk)]';
even_w = 1:2:(exchange.Kk-1); sze_w = size(even_w, 2);
odd_w = 2:2:(exchange.Kk-1); szo_w = size(odd_w, 2);

% for swaps between workers
K_b = params.K - params.W; pair_b = [1:(K_b-1);  2:(K_b)]';
even_b = 1:2:(K_b-1); sze_b = size(even_b, 2);
twr_even = repmat(1:ceil(K_b/2), 2, 1); twr_even = twr_even(:);
odd_b = 2:2:(K_b-1); szo_b = size(odd_b, 2);
twr_odd = repmat(1:ceil(K_b/2), 2, 1); twr_odd=[0; twr_odd(:)];
tww = params.tww; 

% initialise theta and x for each chain
theta = zeros(exchange.Kk, 3, params.W);
x = zeros(exchange.Kk, nET, params.W);
TM = zeros(params.N, params.W+2); tidx=1;
deadline_init.T = 100; deadline.T = params.T;
deadline_init.all = tic; 
for kk=1:exchange.Kk    
    th = params.S_rej(randsample(1:length(params.S_rej), params.W),:); %exprnd(repmat(1./prior, params.W, 1));
    theta(kk,:,:) = th';
    x(kk,:,:) = get_x_multi(params.W, LV, observations, th, deadline_init)';
end
clear('deadline_init', 'idx');

% pre-allocate some arrays
% less tidy for computations
Theta = zeros(exchange.Kk, 3, params.N, params.W); Theta(:, :, 1, :) = theta;
X = zeros(exchange.Kk, nET, params.N, params.W); X(:, :, 1, :) = x;
Rej = zeros(exchange.Kk, params.N, params.W);

% for exchange moves
exK = zeros(exchange.Kk-1, params.W); thetaE = zeros(exchange.Kk-1, 3, params.W);
xE = zeros(exchange.Kk-1, nET, params.W); nE = zeros(exchange.Kk-1, params.W);
epsilonE = zeros(exchange.Kk-1, params.W);

% others
n = zeros(exchange.Kk, params.W); ne = n; nsw = 0; sw = 0; 
resume = zeros(1, params.W);
theta_pending = zeros(2, 3, params.W);
x_pending = zeros(2, 10, params.W);
worK = ones(1, params.W); resumesim = [];
for w=1:params.W
    resumesim(w).resume = 0; %#ok<*AGROW>
    resumesim(w).star = [];
    resumesim(w).xmat = [];
    resumesim(w).target = [];
    resumesim(w).i = []; % simulation will be resumed
    resumesim(w).x = [];
    resumesim(w).tt= [];
    resumesim(w).evolve = [];    
end
burning = 1; evenodd = 1;
within = ones(1, params.W); 
rescount = zeros(exchange.Kk, params.W);

%start the clock
adjust_time = 125/108*1e-5; % to adjust for the different way of measuring time
deadline.all = now; % time checkpoint
deadline.burnin = deadline.all + (params.burnin*adjust_time);
deadline.end = deadline.all + (params.T*adjust_time); 

while now<=deadline.end
    t_c = zeros(1e4, params.W); identifier = zeros(1e4, params.W);
    % next time target on the clock
    deadline.target = min(now + (exchange.deltat)*adjust_time, deadline.end);
    
    %% PARALLEL MOVES %%
    exchange.pair_w = pair_w; exchange.evenodd = evenodd;
    a = tic;
    parfor w=1:params.W
        % Initialise worker
        deadline1 = deadline; tidx1 = 1; id1 = identifier(:,w);
        exchange1 = exchange; rescount1 = rescount(:, w);
        theta1 = theta(:,:,w); x1 = x(:,:,w); kk = worK(w);
        theta_pending1 = theta_pending(:,:,w); x_pending1 = x_pending(:,:,w);
        resumesim1 = resumesim(w);
        within1 = within(w); t_c1 = t_c(:, w);
        
        Theta1 = Theta(:, :, :, w); X1 = X(:, :, :, w);
        Rej1 = Rej(:, :, w); n1 = n(:, w);
        epsilon1 = epsilon(:,w); SIGMA1 = SIGMA(:,:,w);
        params1 = params; resume1 = resume(w);
        % fprintf('\n Worker %d: ', w)
        while(now<=deadline1.target)
            
            %% LOCAL MOVES %%
            % alternate Kk local moves and 1 local exchange move (as per stanndard algorithm)
            lmoves = 0; a1=tic; 
            while (now <= min(deadline1.target, deadline1.end))&&lmoves<exchange.Kk 
                if(~resume1)
                    params1.theta_current = theta1(kk,:);
                    params1.x_current = x1(kk,:);
                else
                    params1.theta_current = theta_pending1;
                    params1.x_current = x_pending1;
                end
                %  fprintf('Worker %d: epsilon = %.3f, SIGMA = ', w, epsilon1(kk))
                
                [params1.theta_current, params1.x_current, resume1, resumesim1, rej1] = LV_local_moves_clock(deadline1, observations, params1, LV, epsilon1(kk), SIGMA1(:,kk), resume1, resumesim1);
                rescount1(kk) = rescount1(kk) + resume1;
                
                if(~resume1)
                    n1(kk) = n1(kk) + 1;
                    theta1(kk,:) = params1.theta_current;
                    x1(kk,:) = params1.x_current;
                    Theta1(kk,:, n1(kk)) = params1.theta_current;
                    X1(kk, :, n1(kk)) = params1.x_current;
                    Rej1(kk, n1(kk)) = rej1;
                    rescount1(kk) = 0; % reset count
                    lmoves = lmoves + (rej1~=2); % one more local move (that wasn't immediately rejected)                    
                    % switch to next chain on worker w
                    kk = kk + 1;
                    if(kk>exchange.Kk)
                        kk = 1;
                    end
                elseif rescount1(kk)>0 %if same chain resumes more than 0 times
                    resume1 = 0; resumesim1.resume=0; % try new proposal
                    rescount1(kk) = 0; % reset count                    
                else
                    theta_pending1 = params1.theta_current;
                    x_pending1 = params1.x_current;                    
                end
            end
            t_c1(tidx1) = toc(a1); id1(tidx1) = 1;              
            tidx1 = tidx1+1;
                     
            
            %% WITHIN WORKER EXCHANGE MOVES %%
            if(now <= min(deadline1.target, deadline1.end))
                b1=tic; 
                
                % exchange move here is standard (does not need to correct
                % for bias)
                if exchange1.evenodd
                    exchange1.evenodd=0;
                    swap_w = exchange1.pair_w(even_w, :);
                    sz = sze_w;
                else
                    exchange1.evenodd=1;
                    swap_w = exchange1.pair_w(odd_w, :);
                    sz = szo_w;
                end
                [theta1, x1, rej_swap, n1, swap_w] = LV_exchange_moves_m(theta1, x1, n1, observations, epsilon1', sz, swap_w);
                % fprintf('| swap chains %d and %d ', swap_w(1,:), swap_w(2,:));
                % if(rej_swap==3)
                %   fprintf('...failed');
                % end
                
                % record swap(s)
                for i=1:sz
                    Rej1(swap_w(i, 1), n1(swap_w(i, 1))) = rej_swap(i);
                    Rej1(swap_w(i, 2), n1(swap_w(i, 2))) = rej_swap(i);
                    Theta1(swap_w(i, 1), :, n1(swap_w(i, 1))) = theta1(swap_w(i, 1),:);
                    Theta1(swap_w(i, 2), :, n1(swap_w(i, 2))) = theta1(swap_w(i, 2),:);
                    X1(swap_w(i, 1), :, n1(swap_w(i, 1))) = x1(swap_w(i, 1), :);
                    X1(swap_w(i, 2), :, n1(swap_w(i, 2))) = x1(swap_w(i, 2), :);
                end
                t_c1(tidx1) = toc(b1); 
                id1(tidx1) = 2; tidx1 = tidx1+1;                
            end
        end
        
        % for between worker exchange moves,
        exK1 = [1:(kk-1) (kk+1):exchange.Kk]; exK(:, w) = exK1;        
        thetaE(:, :, w) = theta1(exK1,:);
        xE(:, :, w)= x1(exK1,:);
        nE(:, w) =  n1(exK1);
        epsilonE(:, w) = epsilon1(exK1);
        
        % updating for every sub-chain kk
        Theta(:,:,:,w) = Theta1; theta(:, :, w) = theta1;
        X(:,:,:,w) = X1; x(:, :, w) = x1;
        Rej(:,:,w) = Rej1; n(:,w) = n1;     
        
        %for resuming after exchange moves
        worK(w) = kk; identifier(:, w) = id1;       
        resume(w) = resume1; resumesim(w) = resumesim1;        
        theta_pending(:,:,w) = theta_pending1;
        x_pending(:,:,w) = x_pending1; 
        within(w) = within1; rescount(:, w) = rescount1;
        t_c(:, w) = t_c1;      
    end
    t_a = toc(a); [sz, max_loc] = max(sum(t_c>0)); tidx_new = tidx+sz-1;
    TM(tidx:tidx_new,:) = [repmat(t_a, [sz, 1]) t_c(1:sz,:) identifier(1:sz,max_loc)];
    tidx = tidx_new + 1; within = ones(1, params.W);
    
    %%%%%%%%%%%%%%%%%%%%%%
    % keeping track of time and burn-in
    if burning&&(now>deadline.burnin)
        ne = n;
        burning=0;
    end    
    TimeRemaining = etime(datevec(deadline.end), clock);
    msg = sprintf('Time remaining: %3.1f seconds', max(0, TimeRemaining)); 
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    %%%%%%%%%%%%%%%%%%%%%%
    
    %% BETWEEN WORKER EXCHANGE MOVES %%        
    if(now<deadline.end)        
        % perform exchange moves between workers before resuming update        
        % select pairs to attempt to swap
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
        theta_swap = [reshape(thetaE(:, 1, :), K_b, 1), reshape(thetaE(:, 2, :), K_b, 1), reshape(thetaE(:, 3, :), K_b, 1)];
        x_swap = zeros(K_b, nET);
        for i=1:nET
            x_swap(:,i) = reshape(xE(:, i, :), K_b, 1);
        end
        epsilon_swap = epsilonE(:)';
        n_swap = nE(:);
        
        % perform exchange moves
        b = tic;
        [theta_swap, x_swap, rej_swap, n_swap, ~] = LV_exchange_moves_m(theta_swap, x_swap, n_swap, observations, epsilon_swap, sz, swap_b);
        TM(tidx,:) = [toc(b) zeros(1, params.W) 3]; tidx = tidx + 1;        
        sw = sw + sum(rej_swap==4);
        
        % redistribute & record
        nE = reshape(n_swap, exchange.Kk-1, params.W);
        exK_ = exK(:);
        for i=swap_b(:)'
            % redistribute
            theta(exK_(i), :, tww(i)) = theta_swap(i, :);
            x(exK_(i), :, tww(i)) = x_swap(i,:);
            n(exK_(i), tww(i)) = n_swap(i);
        
            % record
            Rej(exK_(i), n_swap(i), tww(i)) = rej_swap(twr(i));
            Theta(exK_(i), :, n_swap(i), tww(i)) = theta_swap(i,:);
            X(exK_(i), :, n_swap(i), tww(i)) = x_swap(i,:);
        end                
    end 
end

% tidy up arrays for output
Theta_out = zeros(params.K, 3, params.N);
Rej_out = zeros(params.K, params.N);
X_out = zeros(params.K, nET, params.N);
k = 1;
for w=1:params.W
    for kk=1:exchange.Kk
        X_out(k,:, :) = X(kk, :, :, w);
        Theta_out(k, :, :) = Theta(kk, :, :, w);
        Rej_out(k, :) = Rej(kk, :, w);
        k = k+1;
    end
end

n_out = n(:); ne_out = ne(:);
nmax = max(n_out);
Theta_out = Theta_out(:,:, 1:nmax);
X_out = X_out(:,:, 1:nmax);
Rej_out = Rej_out(:, 1:nmax);
TM = TM(TM(:,end)~=0, :);
end
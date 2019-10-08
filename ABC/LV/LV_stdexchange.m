%% Initialisation
% % Applying approximate Bayesian computation (ABC) with parallel tempering
% % Monte Carlo (PTMC) to % estimate the parameters $\theta =
% % \left(\theta_1, \theta_2, % \theta_3\right)$ of a stochastic
% % Lotka-Volterra predator-prey model.
%
% % Various flags to decide which algorithms are run
% run_single_stdexchange=1; run_multi_stdexchange=1;
% run_single_stdexchange_on_multi=0; run_multi_stdexchange_twice=0;
% true_theta = [1 0.005 0.6];
%
% % observations
% observations.y = [88 165 274 268 114 46 32 36 53 92];
% observations.ET=10+1; observations.dt=1;
%
% % Lotka-Volterra model settings
% LV.M=[50; 100];
% LV.Pre = [1 0; 1 1; 0 1];
% LV.Post = [2 0; 0 2; 0 0];
% LV.h = @(y, th) [th(1)*y(1), th(2)*y(1)*y(2), th(3)*y(2)];
%
% % Run time of algorithm
% params.burnin=3600;
% params.T=params.burnin+4*3600; params.N=1e5; params.K=K;
% fprintf('I will be done on:')
% disp(datetime(clock + [0 0 0 0 0 (run_multi+run_anytime+2)*params.T]))
% % prior
% params.prior = [1 1 1]; %[1 100 1]; %[1 0.01 1];
% params.pr = @(th) params.prior(2)*exp(-params.prior*th');


%% Single processor standard ABC-MCMC algorithm with exchange moves (ABC-PTMC-1)
% The above algorithm but run sequentially on $K=8$ chains, each targetting a posterior
% associated with a different $\varepsilon^k$. The first, _cold_ chain is the
% one associated with ball of radius $\varepsilon^1 = 1$, the _warmer_
% chains target posteriors associated with balls of larger radii. Every
% $\delta_t$ local moves an exchange move occurs between two uniformly
% randomly chosen adjacent chains.

% if(run_single_stdexchange)
%     % Set random seed
%     rng(103505)
%     
%     %     % Number K of chains and corresponding ball radii (epsilon)
% %     K=8; params.K=K;
% %     eps1 = logspace(log10(1), log10(2), fix(K/2));
% %     eps2 = linspace(eps1(end), 10, K+1-fix(K/2)); %logspace(log10(1), log10(10), K);
% %     params.epsilon = [eps1 eps2(2:end)];
%     nocold=0; % allow local moves to occur on the cold chain
%     
%     fprintf('Standard ABC algorithm with %d chains and exchange moves on a single processor\n', params.K)
%     
%     % Exchange move parameters
%     exchange.ipairs = 1:(params.K-1);
%     exchange.pair = [1:(params.K-1);  2:params.K]';
%     exchange.deltat = fix(params.K); % exchange moves occur every K local moves
%     fprintf('Exchange moves every %d local moves for %d seconds\n', exchange.deltat, params.T)
%     
%     % Start from posterior (using rejection ABC as reference
%     params.theta_in = S_rej(randsample(1:length(S_rej), params.K),:); %exprnd(repmat(1./prior, K, 1));
%     
%     
%     % Proposal distribution covariance varies from chain to chain
% %     sigma = [0.01 0.02 0.05 0.1 0.2 0.3 0.4 0.5]; %[0.0015 0.004 0.015 0.04 0.075 0.15 0.3 0.5] %[0.004 0.015 0.04 0.066 0.075 0.15 0.18 0.25];
% %     params.SIGMA = [sigma; sigma*1e-2; sigma];
%     
%     % Run ABC-PTMC-1 algorithm
%     tic
%     [Theta_e, X_e, ne_e, n_e, Rej_e, nsw, sw, TM] = LV_ABC_standard_exchange(params, exchange, LV, observations, nocold);
%     T_e = toc %#ok<*NOPTS>
%     
%     % Record cold chain
%     k=1;
%     S_e = Theta_e(k, :, ne_e(k):n_e(k)); S_e = reshape(S_e, size(S_e, 2), size(S_e,3))';
% %     dlmwrite(sprintf('results/ABC/LV/LV_standard_exchange_%d_%d.csv', n_e(k)-ne_e(k)+1, params.T), S_e);
%     %  file = sprintf('results/ABC/LV/LV_standard_exchange_%d_%d_X.mat', n_e(k)-b_e(k)+1, params.T);
%     %  dlmwrite(file, X_e);
%     
%     % Break down the composition of each chain and compute local and
%     % exchange moves acceptance rates (as percentages)
%     [out_e, rates_e, b_e] = print_summary_chain({Rej_e, n_e, ne_e, nsw, sw}, {params.epsilon, params.SIGMA}, 0, 1);
% 
%    
% end

%% Single processor standard ABC-MCMC algorithm with exchange moves (ABC-PTMC-1) W times in parallel
% The above algorithm but run W times in parallel on W processors, i.e.
% equivalent to W separate runs of the algorithm.

if(run_single_stdexchange_on_multi)
    % Set random seed
    rng(103505)
    
    % Start parallel pool on 4 processors
%     W = 4;
    if(isempty(gcp('nocreate')))
        p = parpool(min(W, 4));
        p.IdleTimeout = 9000;
    end
    
    % Number K of chains and corresponding ball radii (epsilon)
    %K=6; params.K=K;
    %eps1 = logspace(log10(1), log10(1.5), fix(params.K-2));
    %eps2 = [1.5 11 15];%linspace(eps1(end), 15, K+1-fix(K-2)); %logspace(log10(1), log10(10), K);
    %params.epsilon = [eps1 eps2(2:end)];
    nocold=0; % allow local moves to occur on the cold chain
    
    fprintf('Single processor standard ABC algorithm with %d chains and exchange moves \n', params.K)
    fprintf('---> run %d times in parallel \n', W)
    
    % Exchange move parameters
    exchange.ipairs = 1:(params.K-1);
    exchange.even = 1:2:(params.K-1);
    exchange.odd = 2:2:(params.K-1);
    exchange.pair = [1:(params.K-1);  2:params.K]';
    exchange.deltat = fix(params.K*nsweeps); % exchange moves occur every K local moves
    fprintf('Exchange moves every %d local moves for %d seconds\n', exchange.deltat, params.T) 
    
    
    % Proposal distribution covariance varies from chain to chain
    %sigma = [0.008 0.025 0.05 0.09 0.25 0.5];% [0.008 0.025 0.05 0.09 0.25 0.4];
    params.SIGMA = [sigma; sigma*1e-2; sigma];
    
    % For parfor loop
    n_e = zeros(params.K, W); ne_e =n_e;
    epsilon_out_e = zeros(params.K, W);
    nsw_e = zeros(1, W); sw_e = nsw_e;
    Theta_e = cell(1, W); %zeros(K, 3, params.N, W);
    X_e = cell(1, W); %zeros(K, observations.ET-1, params.N, W) ;
    Rej_e = cell(1, W); %zeros(K, params.N, W) ;
    TM = cell(1, W); median_times = zeros(1, W); %for ABC-APTMC
    S_e = cell(1, W);
    theta_in_e = cell(1, W);
    % Start from posterior (using rejection ABC as reference)
    for w=1:W
        theta_in_e{w} = S_rej(randsample(1:length(S_rej), params.K),:); %exprnd(repmat(1./prior, K, 1));
    end
    
    % Run ABC-PTMC-1 algorithm
    tic
    parfor w=1:W
        params1 = params;
        params1.theta_in = theta_in_e{w};
        
        [Theta_e1, X_e1, ne_e1, n_e1, Rej_e1, nsw1, sw1, TM1] = LV_ABC_standard_exchange(params1, exchange, LV, observations, nocold);
    
        Theta_e{w} = Theta_e1;
        X_e{w} = X_e1;
        Rej_e{w} = Rej_e1;
        TM{w} = TM1; median_times(w) = median(TM1(TM1(:,2)==1,1));
        nsw_e(w) = nsw1; sw_e(w) = sw1;
        n_e(:, w) = n_e1; ne_e(:,w) = ne_e1;
        %epsilon_out_e(:,w) = epsilon1;
        
        % Isolate cold chain(s)
        k1=1; S_e1 = Theta_e1(k1, :, ne_e1(k1):n_e1(k1)); 
        S_e1 = reshape(S_e1, size(S_e1, 2), size(S_e1,3))';
        S_e{w} = S_e1;
    end
    T_e = toc %#ok<*NOPTS>
    
    k=1; b_e = zeros(params.K, W); %n_e - ne_e+1; 
    out_e = zeros(params.K, 5, W);
    rates_e = zeros(params.K, 2, W);
    for w=1:W
        
        % Break down the composition of each chain and compute local and
        % exchange moves acceptance rates (as percentages)
        fprintf(' \n --------------- WORKER %d --------------- \n \n', w)
        [out_e(:,:,w), rates_e(:,:,w), b_e(:, w)] = print_summary_chain({Rej_e{w}, n_e(:,w), ne_e(:,w), nsw_e(w), sw_e(w)}, {params.epsilon, params.SIGMA}, 0, 1);
        
        
        % Record cold chain(s)
        %  dlmwrite(sprintf('results/ABC/LV/LV_standard_exchange_%d_%d_%d.csv', b_e(k, w), params.T), S_e{w}, w);
        %  file = sprintf('results/ABC/LV/LV_standard_exchange_%d_%d_X.mat', b_e(k, w), params.T);
        %  dlmwrite(file, X_e{w});
        
    end
    
    
end


%% Multi-processor standard ABC-MCMC algorithm with exchange moves (ABC-PTMC-W)
% The above algorithm but run in parallel on 4 workers, with $Kk=2$
% chains per processor (i.e. K=8 chains in total), each chain targets a
% posterior associated with a different $\varepsilon^k$ and exchange moves
% are performed between workers.

% if(run_multi_stdexchange)
%     % Run on 4 workers
%     params.W = 2; params.Ww=2; exchange.Kk=3;%#ok<*UNRCH>
%     
%     fprintf('Standard ABC algorithm with %d chains, i.e. %d chains per worker on %d workers\n', exchange.Kk*params.W, exchange.Kk, params.W)
%     params.K=params.W*exchange.Kk; % exchange.deltat = fix(5);
%     exchange.deltat = fix(exchange.Kk*2); % exchange moves occur every K local moves
%     fprintf('Exchange moves every %d local moves for %d seconds\n', exchange.deltat, params.T)
%     
%     %Adapt epsilon and sigma to number of chains
%     sigma = [0.01 0.025 0.05 0.1 0.25 0.5]; %[0.008 0.01 0.03 0.05 0.075 0.09 0.15 0.25 0.35 0.45 0.55 0.65];%linspace(0.005, 0.4, params.K); %[0.004 0.015 0.04 0.066 0.075 0.15 0.18 0.25];
%     params.SIGMA = [sigma; sigma*1e-2; sigma];
%     eps1 = logspace(log10(1), log10(2), fix(params.K-1));
%     eps2 = [2 11];%15];%linspace(eps1(end), 15, K+1-fix(K-2)); %logspace(log10(1), log10(10), K);
%     params.epsilon = [eps1 eps2(2:end)];
%     
%     % Start parallel pool
%     if(isempty(gcp('nocreate')))
%         p = parpool(min(W, 4));
%         p.IdleTimeout = 9000;
%     end
%     multi_temp_per_worker = 1; % multiple temperatures on each worker
%     
%     % Run ABC-PTMC-4 algorithm
%     % Set random seed
%     rng(65416)
%     tic
%     [Theta_em, X_em, n_em, ne_em, Rej_em, nswm, swm, TM_m] = LV_ABC_standard_exchange_multi_2(params, exchange, LV, observations, 0);
%     T_em = toc
%     
%     % Record cold chain
%     k=1;
%     S_em = Theta_em(k, :, ne_em(k):n_em(k)); S_em = reshape(S_em, size(S_em, 2), size(S_em,3))';
%     % dlmwrite(sprintf('results/ABC/LV/LV_standard_exchange_multi_%d.csv', n_em(k)-b_em(k)+1), S_em);
%     
%     % Break down the composition of each chain and compute local and
%     % exchange moves acceptance rates (as percentages)
%     [out_em, rates_em, b_em] = print_summary_chain({Rej_em, n_em, ne_em, nswm, swm}, {params.epsilon, params.SIGMA}, 1, 1);
%     out_em(:, 6, :) = sum(out_em(:, [3 5], :), 2);
% 
% end

%% Multi-processor standard ABC-MCMC algorithm with exchange moves (ABC-PTMC-W)
% The above algorithm but run in parallel on 4 workers, with $Kk=2$
% chains per processor (i.e. K=8 chains in total), each chain targets a
% posterior associated with a different $\varepsilon^k$ and exchange moves
% are performed between workers.

if(run_multi_stdexchange_multi)
    % Run on 2 workers twice
    rng(65416)
    params.W = 4; params.Ww=4; exchange.Kk=5;%#ok<*UNRCH>
    
    fprintf('Standard ABC algorithm with %d chains, i.e. %d chains per worker on %d workers\n', exchange.Kk*params.W, exchange.Kk, params.W)
    params.K=params.Ww*exchange.Kk; % exchange.deltat = fix(5);
    exchange.deltat = fix(exchange.Kk*nsweeps); % exchange moves occur every K local moves
    fprintf('Exchange moves every %d local moves for %d seconds\n', exchange.deltat, params.T)
    
    %Adapt epsilon and sigma to number of chains
    s1 = logspace(log10(0.008), log10(0.06), fix(params.K-exchange.Kk));%[0.008 0.015 0.025 0.05 0.075 0.1 0.5 1];%[0.01 0.05 0.15 0.2 0.5 1]; %[0.008 0.025 0.05 0.09 0.25 0.5];
    s2 = logspace(log10(s1(end)), log10(.5), fix(exchange.Kk)+1);
%     s3 = linspace(0.08, .5, fix(exchange.Kk));
    sigma = [s1 s2(2:end)];
    params.SIGMA = [sigma; sigma*1e-2; sigma];
    eps1 = logspace(log10(1), log10(1.5), fix(params.K-2*exchange.Kk));
    eps2 = logspace(log10(eps1(end)), log10(2.5), fix(exchange.Kk)+1);
    eps3 = logspace(log10(eps2(end)), log10(15), fix(exchange.Kk)+1); %logspace(log10(1), log10(10), K);
    params.epsilon = [eps1 eps2(2:end) eps3(2:end)];
    params.tww = repmat(1:params.W, exchange.Kk, 1);params.tww=params.tww(:);
    params.twk = repmat(1:exchange.Kk, 1, params.W);
    params.twr = repmat(1:ceil(params.K/2), 2, 1); params.twr=params.twr(:);
    
    % Start parallel pool
    if(isempty(gcp('nocreate')))
        p = parpool(min(params.W, 4));
        p.IdleTimeout = 9000;
    end
    multi_temp_per_worker = 1; % multiple temperatures on each worker
    
    % Run ABC-PTMC-4 algorithm
    % Set random seed
     % For parfor loop
    n_em = zeros(params.K, ntimes); ne_em = n_em;
    nswm = zeros(1, ntimes); swm = nswm;
    Theta_em = cell(1, ntimes); X_em = cell(1, ntimes); 
    Rej_em = cell(1, ntimes); median_times = zeros(1, ntimes);
    TM_m = cell(1, ntimes); S_em = cell(1, ntimes); 
    TMall_m = cell(1,ntimes);
    
    % To break down the composition of each chain and compute local and
    % exchange moves acceptance rates (as percentages)
    b_em = zeros(params.K, ntimes); %an_e - ane_e+1;
    out_em = zeros(params.K, 5, ntimes);
    rates_em = zeros(params.K, 2, ntimes);
    k=1;
    
    tic
    for idx=1:ntimes
        fprintf(' \n --------------- RUN No. %d --------------- \n \n', idx)
        % Run algorithm
        [Theta_em{idx}, X_em{idx}, n_em(:, idx), ne_em(:, idx), Rej_em{idx}, nswm(idx), swm(idx), TM1] = LV_ABC_standard_exchange_multi_2(params, exchange, LV, observations, 0);
        
        %TM1 = TM_M{idx};
        TMall1 = get_times(TM1);
        median_times(idx) = median(TMall1(:,2));%max(median(TM1(TM1(:,params.W+2)==1, 2:(params.W+1))));
        
        TM_m{idx} = TM1;
        TMall_m{idx} = TMall1;
        % Record cold chain
        S_em{idx} = Theta_em{idx}(k, :, ne_em(k, idx):n_em(k, idx)); S_em{idx} = reshape(S_em{idx}, size(S_em{idx}, 2), size(S_em{idx},3))';
        [out_em(:,:,idx), rates_em(:,:,idx), b_em(:,idx)] = print_summary_chain({Rej_em{idx}, n_em(:, idx), ne_em(:, idx), nswm(idx), swm(idx)}, {params.epsilon, params.SIGMA}, 0, 1);
        
    end
    T_em = toc
    
    out_em(:, 6, :) = sum(out_em(:, [3 5], :), 2);
    clear('TM1', 'TMall1')

end
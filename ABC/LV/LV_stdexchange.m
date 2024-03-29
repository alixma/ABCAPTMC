%% Single processor standard ABC-MCMC algorithm with exchange moves (ABC-PTMC-1) W times in parallel
% ABC-MCMC algorithm un sequentially on $K$ chains, each targetting a posterior
% associated with a different $\varepsilon^k$. The first, _cold_ chain is the
% one associated with ball of radius $\varepsilon^1 = 1$, the _warmer_
% chains target posteriors associated with balls of larger radii. Every
% $\delta_t$ local moves an exchange move occurs between two uniformly
% randomly chosen adjacent chains. It ist run W times in parallel on W 
% processors, i.e. equivalent to W separate runs of the algorithm.

if(run_single_stdexchange)
    % Set random seed
    rng(103505)
    
    if ~run_multi % do not initialise in multi processor example
        W = 4;  
        % Start parallel pool on 4 processors
        if(isempty(gcp('nocreate')))
            p = parpool(min(W, 4));
            p.IdleTimeout = 9000;
        end
        
        % Number K of chains and corresponding ball radii (epsilon)        
        K = 6; params.K = K;
        eps1 = logspace(log10(1), log10(1.5), fix(params.K-2));
        eps2 = [1.5 11 15];%linspace(eps1(end), 15, K+1-fix(K-2)); %logspace(log10(1), log10(10), K);
        params.epsilon = [eps1 eps2(2:end)];
        clear('eps1', 'eps2');
        
        % Proposal distribution covariance varies from chain to chain
        sigma = [0.008 0.025 0.05 0.09 0.25 0.5];
        params.SIGMA = [sigma; sigma*1e-2; sigma];
    end
    
    % Exchange move parameters
    exchange.pair = [1:(params.K-1);  2:(params.K)]';
    exchange.deltat = fix(params.K*nsweeps); % exchange moves occur every K*nsweeps local moves
    
    fprintf('Single processor standard ABC algorithm with %d chains and exchange moves \n', params.K)
    fprintf('---> run %d times in parallel \n', W)        
    fprintf('Exchange moves every %d local moves for %d seconds\n', exchange.deltat, params.T)
    
    % For parfor loop
    n_e = zeros(params.K, W); ne_e = n_e;
    nsw_e = zeros(1, W); sw_e = nsw_e;
    Theta_e = cell(1, W);  X_e = cell(1, W);
    Rej_e = cell(1, W); S_e = cell(1, W);
    TM = cell(1, W); median_times = zeros(1, W); %for ABC-APTMC-1    
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
        
        [Theta_e1, X_e1, Rej_e1, n_e1, ne_e1, nsw1, sw1, TM1] = LV_ABC_standard_exchange_1(params1, exchange, LV, observations);
    
        Theta_e{w} = Theta_e1; X_e{w} = X_e1; Rej_e{w} = Rej_e1;      
        TM{w} = TM1; median_times(w) = median(TM1(TM1(:,2)==1,1));
        nsw_e(w) = nsw1; sw_e(w) = sw1;
        n_e(:, w) = n_e1; ne_e(:,w) = ne_e1;
        
        % Isolate cold chain(s)
        k1 = 1; S_e1 = Theta_e1(k1, :, ne_e1(k1):n_e1(k1)); 
        S_e1 = reshape(S_e1, size(S_e1, 2), size(S_e1,3))';
        S_e{w} = S_e1;
    end
    T_e = toc %#ok<*NOPTS>
    
    k=1; b_e = zeros(params.K, W); 
    out_e = zeros(params.K, 6, W); rates_e = zeros(params.K, 2, W);    
    for w=1:W        
        % Break down the composition of each chain and compute local and
        % exchange moves acceptance rates (as percentages)
        fprintf(' \n --------------- WORKER %d --------------- \n \n', w)
        [out_e(:,:,w), rates_e(:,:,w), b_e(:, w)] = print_summary_chain({Rej_e{w}, n_e(:,w), ne_e(:,w), nsw_e(w), sw_e(w)}, {params.epsilon, params.SIGMA}, 0, 1);
        
        % Record cold chain(s)
        %  dlmwrite(sprintf('results/LV/LV_standard_exchange_1_%d_%d_%d.csv', b_e(k, w), params.T), S_e{w}, w);
        %  file = sprintf('results/LV/LV_standard_exchange_1_%d_%d_X.mat', b_e(k, w), params.T);
        %  dlmwrite(file, X_e{w});        
    end    
    
end

%% Multi-processor standard ABC-MCMC algorithm with exchange moves (ABC-PTMC-W)
% The above algorithm but run in parallel on W workers, with $Kk=2$
% chains per processor (i.e. K=8 chains in total), each chain targets a
% posterior associated with a different $\varepsilon^k$ and exchange moves
% are performed between workers.

if(run_multi_stdexchange)
    rng(65416)
    params.W = 4; exchange.Kk = 5;
    
    fprintf('Standard ABC algorithm with %d chains, i.e. %d chains per worker on %d workers\n', exchange.Kk*params.W, exchange.Kk, params.W)
    params.K = params.W*exchange.Kk; % exchange.deltat = fix(5);
    exchange.deltat = fix(exchange.Kk*nsweeps); % exchange moves occur every K local moves
    fprintf('Exchange moves every %d local moves for %d seconds\n', exchange.deltat, params.T)
    
    % Adapt epsilon and sigma to number of chains
    s1 = logspace(log10(0.008), log10(0.06), fix(params.K-exchange.Kk));
    s2 = logspace(log10(s1(end)), log10(.5), fix(exchange.Kk)+1);
    sigma = [s1 s2(2:end)];
    params.SIGMA = [sigma; sigma*1e-2; sigma];
    eps1 = logspace(log10(1), log10(1.5), fix(params.K-2*exchange.Kk));
    eps2 = logspace(log10(eps1(end)), log10(2.5), fix(exchange.Kk)+1);
    eps3 = logspace(log10(eps2(end)), log10(15), fix(exchange.Kk)+1); 
    params.epsilon = [eps1 eps2(2:end) eps3(2:end)];
    clear('s1', 's2', 'eps1', 'eps2', 'eps3')
    
    % various exchange move parameters
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
    b_em = zeros(params.K, ntimes); 
    out_em = zeros(params.K, 6, ntimes);
    rates_em = zeros(params.K, 2, ntimes);
    k=1;
    
    tic
    for idx=1:ntimes
        fprintf(' \n --------------- RUN No. %d --------------- \n \n', idx)
        % Run algorithm
        [Theta_em{idx}, X_em{idx}, n_em(:, idx), ne_em(:, idx), Rej_em{idx}, nswm(idx), swm(idx), TM1] = LV_ABC_standard_exchange_W(params, exchange, LV, observations);
        
        %TM1 = TM_M{idx};
        TMall1 = get_times(TM1);
        median_times(idx) = median(TMall1(:,2));
        
        TM_m{idx} = TM1;
        TMall_m{idx} = TMall1;
        % Record cold chain
        S_em{idx} = Theta_em{idx}(k, :, ne_em(k, idx):n_em(k, idx)); S_em{idx} = reshape(S_em{idx}, size(S_em{idx}, 2), size(S_em{idx},3))';
        [out_em(:,:,idx), rates_em(:,:,idx), b_em(:,idx)] = print_summary_chain({Rej_em{idx}, n_em(:, idx), ne_em(:, idx), nswm(idx), swm(idx)}, {params.epsilon, params.SIGMA}, 0, 1);
        
    end
    T_em = toc
    clear('TM1', 'TMall1')

end
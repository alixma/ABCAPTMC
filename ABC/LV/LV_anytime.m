%% Initialisation
% % Applying anytime parallel tempering Monte Carlo (APTMC) to
% % estimate the parameters $\theta = \left(\theta_1, \theta_2,
% % \theta_3\right)$ of a stochastic Lotka-Volterra predator-prey model.
% 
%
% % Various flags to decide which algorithms are run
% run_single_anytime=1; run_multi_anytime=1;
% run_single_anytime_on_multi=0;
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
% disp(datetime(clock + [0 0 0 0 0 (run_single_anytime+run_multi_anytime)*params.T]))
% % prior
% params.prior = [1 1 1]; %[1 100 1]; %[1 0.01 1];
% params.pr = @(th) params.prior(2)*exp(-params.prior*th');


%% Single processor anytime ABC-MCMC algorithm with exchange moves (ABC-APTMC-1)
% The anytime version of the ABC_APTMC-1 algorithm, corrected for length bias.
% Again run on $K=8$ chains, each targetting a posterior associated with a
% different $\varepsilon^k$. This time, exchange moves occur every
% $\delta^a_t$ seconds, where now $\delta^a_t$ is the median time the
% previous algorithm took to perform $\delta_t$ local moves.


% if(run_single_anytime)
%     % Set random seed
%     rng(65786)
%     
%     % Number K of chains and corresponding ball radii (epsilon)
%     % K=8; params.K=K;
%     % eps1 = logspace(log10(1), log10(2), fix(K/2)); 
%     % eps2 = linspace(eps1(end), 10, K+1-fix(K/2)) ; %logspace(log10(1), log10(10), K);
%     % params.epsilon = [eps1 eps2(2:end)];    
%     % nocold=0; % allow local moves to occur on the cold chain
%     
%     % More exchange moves parameters
%     exchange.ipairs = 1:(params.K-2);
%     exchange.pair = [1:(params.K-2);  2:(params.K-1)]';
%     % Ensure median time spent performing local moves before each exchange
%     % move remains the same as for ABC-PTMC
%     exchange.deltat = round(median(TM(TM(:,2)==1,1)),3);
%     
%     fprintf('Anytime ABC algorithm with %d chains and exchange moves on a single processor\n', params.K)    
%     fprintf('Exchange moves every %.2f seconds for %d seconds\n', exchange.deltat, params.T)
%            
%     % Start from same place as ABC-PTMC-1
%     params.theta_in = Theta_e(:,:,1); % params.S_rej(randsample(1:length(params.S_rej), K),:); %exprnd(repmat(1./prior, K, 1));
%    
%     
%     %Run ABC-APTMC-1 algorithm
%     tic
%     [Theta_ae, X_ae, Rej_ae, an_e, ane_e, answ, asw, TM_a, alln] = LV_ABC_anytime_all(params.K, observations, LV, params, exchange, nocold);
%     toc
%     
%     % Record cold chain
%     k=1; 
%     A_e = Theta_ae(k,:, ane_e(k):an_e(k)); A_e = reshape(A_e, size(A_e, 2), size(A_e,3))';
%     dlmwrite(sprintf('results/ABC/LV/LV_anytimem_%d_%d.csv', an_e(k)-ane_e(k)+1, params.T), A_e);
%     
%     % Break down the composition of each chain and compute local and
%     % exchange moves acceptance rates (as percentages)
%     [o1, r1, b_a] = print_summary_chain({Rej_ae, an_e, ane_e, answ, asw}, {params.epsilon, params.SIGMA}, 0, 1);
%     
%     alln = alln(alln(:,1)~=0,:);
% end

%% Single processor anytime ABC-MCMC algorithm with exchange moves (ABC-APTMC-1) W times in parallel
% The above algorithm but run W times in parallel on W processors, i.e.
% equivalent to W separate runs of the algorithm.

if(run_single_anytime_on_multi)
    % Set random seed
    rng(65786)
    
    %     % Number K of chains and corresponding ball radii (epsilon)
    % K=8; params.K=K;
    % eps1 = logspace(log10(1), log10(2), fix(K/2));
    % eps2 = linspace(eps1(end), 10, K+1-fix(K/2)) ; %logspace(log10(1), log10(10), K);
    % params.epsilon = [eps1 eps2(2:end)];
    % nocold=0; % allow local moves to occur on the cold chain
    
    %     % Start parallel pool on 4 processors
    %     W = 4;
    %     if(isempty(gcp('nocreate')))
    %         p = parpool(min(W, 4));
    %         p.IdleTimeout = 9000;
    %     end
    
    % More exchange moves parameters
    exchange.ipairs = 1:(params.K-2);
    exchange.pair = [1:(params.K-2);  2:(params.K-1)]';
    % Ensure median time spent performing local moves before each exchange
    % move remains the same as for ABC-PTMC
    exchange.deltat = round(median(median_times),3);% to account for communication overhead
    %params.epsilon;
    %sigma ;
    
    fprintf('Anytime ABC algorithm with %d chains and exchange moves on a single processor\n', params.K)
    fprintf('---> run %d times in parallel \n', W)    
    fprintf('Exchange moves every %.2f seconds for %d seconds\n', exchange.deltat, params.T)
    
    % For parfor loop
    an_e = zeros(params.K, W); ane_e =an_e;
    answ = zeros(1, W); asw = answ;
    epsilon_out_a = zeros(params.K, W);
    Theta_ae = cell(1, W); %zeros(K, 3, params.N, W);
    X_ae = cell(1, W); %zeros(K, observations.ET-1, params.N, W) ;
    Rej_ae = cell(1, W); %zeros(K, params.N, W) ;
    TM_a = cell(1, W); A_e = cell(1, W);
    %     theta_in_a = cell(1, W);
    %     % Start from posterior (using rejection ABC as reference)
    %     for w=1:W
    %         theta_in_a{w} = S_rej(randsample(1:length(S_rej), K),:); %exprnd(repmat(1./prior, K, 1));
    %     end
    
    %Run ABC-APTMC-1 algorithm
    tic
   parfor w=1:W
        params1 = params;
        % Start from same place as ABC-PTMC-1
        params1.theta_in = theta_in_e{w}; %theta_in_a{w};
        
        [Theta1, X1, Rej1, n1, ne1, nsw1, sw1, TM1] = LV_ABC_anytime_all(params.K, observations, LV, params1, exchange, nocold);
        
        Theta_ae{w} = Theta1;
        X_ae{w} = X1;
        Rej_ae{w} = Rej1;
        TM_a{w} = TM1;
        answ(w) = nsw1; asw(w) = sw1;
        an_e(:, w) = n1; ane_e(:,w) = ne1;
        %epsilon_out_a(:,w) = epsilon1;
        
        % Isolate cold chain(s)
        k1=1; A1 = Theta1(k1, :, ne1(k1):n1(k1));
        A1 = reshape(A1, size(A1, 2), size(A1,3))';
        A_e{w} = A1;
    end
    T_a = toc      
    
    k=1; b_a = zeros(params.K, W); %an_e - ane_e+1; 
    out_a = zeros(params.K, 5, W);
    rates_a = zeros(params.K, 2, W);
    for w=1:W
        % Break down the composition of each chain and compute local and
        % exchange moves acceptance rates (as percentages)
        fprintf(' \n --------------- WORKER %d --------------- \n \n', w)
        [out_a(:,:,w), rates_a(:,:,w), b_a(:,w)] = print_summary_chain({Rej_ae{w}, an_e(:,w), ane_e(:,w), answ(w), asw(w)}, {params.epsilon, params.SIGMA}, 0, 1);
        
        % Record cold chain(s)        
        %  dlmwrite(sprintf('results/ABC/LV/LV_anytime_%d_%d_%d.csv', b_a(k, w), params.T), A_e{w}, w);
        %  file = sprintf('results/ABC/LV/LV_anytime_%d_%d_X.mat', b_a(k, w), params.T);
        %  dlmwrite(file, X_ae{w});                        
    end         
    
end

%% Multi-processor anytime ABC-MCMC algorithm with exchange moves (ABC-APTMC-W)
% The above algorithm but run in parallel on 4 workers, with $Kk=2$
% chains per processor (i.e. K=8 chains in total), each chain targets a
% posterior associated with a different $\varepsilon^k$ and exchange moves
% are performed between workers.

if(run_multi_anytime)
    % Set random seed
    rng(405310)    
    
    % Run on 4 workers
    % params.W = 4; exchange.Kk=3;%#ok<*UNRCH>
    
    fprintf('Anytime ABC algorithm with %d chains, i.e. %d chains per worker on %d workers\n', exchange.Kk*params.W, exchange.Kk, params.W)
    exchange.deltat = round(max(median(TM_m(TM_m(:,params.W+2)==1, 2:(params.W+1)))),3); %round(median(TM_m(TM_m(:,params.W+2)==1,1)),3); %2;
    fprintf('Exchange moves every %.2f seconds for around %d seconds\n', exchange.deltat, params.T)
    
    if(isempty(gcp('nocreate')))
        parpool(min(params.W, 4));
    end
    multi_temp_per_worker = 1;
    
    tic
    [Theta_aem, X_aem, Rej_aem, an_em, ane_em, answ_em, asw_em, TM_am] = LV_ABC_anytime_multiclock_all_2(params, exchange, LV, observations, 0);
    T_am = toc
    k=1;
    % Record cold chain
    if(multi_temp_per_worker)
        A_em = Theta_aem(k,:,ane_em(k):an_em(k));
        A_em = reshape(A_em, size(A_em, 2), size(A_em, 3))';
    else
        A_em = [reshape(Theta_em(1,:,1:an_em(1,1)), 3, an_em(1,1))' ones(an_em(1,1), 1);
            reshape(Theta_em(2,:,1:an_em(2,1)), 3, an_em(2,1))' 2*ones(an_em(2,1), 1)];
    end
    
    % Break down the composition of each chain and compute local and
    % exchange moves acceptance rates (as percentages)
    [out_am, rates_am, b_am] = print_summary_chain({Rej_aem, an_em, ane_em, answ_em, asw_em}, {params.epsilon, params.SIGMA}, 1, 1);
    
    out_am(:, 6, :) = sum(out_am(:, [3 5], :), 2);


    
end

%% Multi-processor anytime ABC-MCMC algorithm with exchange moves (ABC-APTMC-W)
% The above algorithm but run twice in parallel on 2 workers, with $Kk=3$
% chains per processor (i.e. K=6 chains in total), each chain targets a
% posterior associated with a different $\varepsilon^k$ and exchange moves
% are performed between workers.

if(run_multi_anytime_multi)
    rng(45310)
    
    % Run on 4 workers
    % params.W = 4; exchange.Kk=3;%#ok<*UNRCH>
    
    fprintf('Anytime ABC algorithm with %d chains, i.e. %d chains per worker on %d workers\n', exchange.Kk*params.W, exchange.Kk, params.W)
    exchange.deltat = round(median(median_times),3); %round(median(TM_m(TM_m(:,params.W+2)==1,1)),3);  %2;
    fprintf('Exchange moves every %.2f seconds for around %d seconds\n', exchange.deltat, params.T)
    
    if(isempty(gcp('nocreate')))
        p = parpool(min(params.W, 4));
        p.IdleTimeout = 9000;
    end
    multi_temp_per_worker = 1;
    
    % For parfor loop
    an_em = ones(params.K, ntimes); ane_em = an_em;
    answ_em = ones(1, ntimes); asw_em = answ_em;
    Theta_aem = cell(1, ntimes); X_aem = cell(1, ntimes); 
    Rej_aem = cell(1, ntimes); 
    TM_am = cell(1, ntimes); A_em = cell(1, ntimes); 
    
    % To break down the composition of each chain and compute local and
    % exchange moves acceptance rates (as percentages)
    b_am = zeros(params.K, ntimes); %an_e - ane_e+1;
    out_am = zeros(params.K, 5, ntimes);
    rates_am = zeros(params.K, 2, ntimes);
    
    k=1;
    tic
    dirname =  sprintf('LV_multi_20181127_%s_%d', 'anytime', params.T);
    %     create folder to store results
    mkdir(sprintf('results/ABC/LV/%s', dirname))
    params.tww = repmat(1:params.W, exchange.Kk-1, 1);params.tww=params.tww(:);
    
    
    for idx=1:ntimes
        fprintf(' \n --------------- RUN No. %d --------------- \n \n', idx)
        
        % Run algorithm
        [Theta_aem{idx}, X_aem{idx}, Rej_aem{idx}, an_em(:, idx), ane_em(:, idx), answ_em(idx), asw_em(idx), TM_am{idx}] = LV_ABC_anytime_multiclock_all_2(params, exchange, LV, observations, 0);
        fprintf('\n')
        
        % Record cold chain
        A_em{idx} = Theta_aem{idx}(k, :, max(1,ane_em(k, idx)):an_em(k, idx)); A_em{idx} = reshape(A_em{idx}, size(A_em{idx}, 2), size(A_em{idx},3))';
        [out_am(:,:,idx), rates_am(:,:,idx), b_am(:,idx)] = print_summary_chain({Rej_aem{idx}, an_em(:, idx), ane_em(:, idx), answ_em(idx), asw_em(idx)}, {params.epsilon, params.SIGMA}, 0, 1);
        
        save_one(idx, dirname, 'anytimem', A_em, b_am, TM_am, out_am, rates_am)
        
    end
    T_am = toc
    
    out_am(:, 6, :) = sum(out_am(:, [3 5], :), 2);
end
%% PERFORMANCE COMPARISON OF ABC-APTMC ALGORITHMS

%% Initialisation

fprintf('Initiating... \n')
% Various parameters for simulating observations of the MA(q) process
simulations.q=2;
theta_true=[0.6; 0.2]; simulations.nsamples=500;
maxNumCompThreads(1);

C=8; 
W = min(C, 4);
if(isempty(gcp('nocreate')))&&(C>1)    
    parpool(W);
end

% Settings of chains for exchange moves
% number of chains for parallel tempering
K=10; params.K = K;
params.correct=1;

% Burnin and total running time of algorithms
b=60; params.T = b+3*3600; params.N=1e6; params.burnin = b;

fprintf('I will be done on:')
disp(datetime(clock + [0 0 0 0 0 (ceil(C/W)*params.T*4)]))

% Define the radius $\varepsilon$ of the ball associated with each chain
eps1 = logspace(log10(0.02), log10(0.1), fix(K/2)+1); 
eps2 = logspace(log10(eps1(end)), log10(1), K+1-length(eps1));
params.epsilon = [eps1, eps2(2:end)];
clear('eps1', 'eps2');

fprintf('Initiating...done \n\n')

%% Observations
% Simulate the observations that will be used throughout this example.

fprintf('Simulating observations... \n')
rng(324721);
% simulate MA process
Y = MA_sim(theta_true, simulations, 1, 0);
% compute the first q sample autocorrelations
[xc, lg] = xcov(Y, simulations.q, 'coef');
simulations.y = xc(lg>0)';
fprintf('y = ')
disp(simulations.y)
clear('xc', 'lg')
fprintf('Simulating observations...done \n\n')

% Posterior from rejection sampling for reference
R_rej = dlmread(sprintf('results/MA/MA_ABC_rejection_%d.csv', 13189));

%% ALGORITHMS runs
fprintf('Running algorithms... \n')
rng('default');

%% Standard ABC-MCMC algorithm
% The vanilla MCMC algorithm which employs ABC in its local moves.
fprintf('Standard ABC for %d seconds... \n', params.T)
params.deltat = params.T+1000; params.K = 1;
% Start from the posterior (obtained by rejection ABC)
params.theta_in = R_rej(randsample(1:length(R_rej), params.K),:)';

% Standard deviation for Gaussian random walk proposal for local moves
params.rho=0.25;

tic
[Theta_s, X_s, Rej_s, n_s, ne_s] = MA_ABC_standard(C, params, simulations, params.epsilon(:,1));
b_s = n_s'-ne_s'+1
toc
[~, c_s] = max(b_s);
fprintf('Standard ABC... done \n\n')

%% Standard ABC-MCMC algorithm with exchange moves (ABC-PTMC)
% The above algorithm but run on $K=8$ chains, each targetting a different posterior.

% Standard deviation for Gaussian random walk proposal for local moves
params.rho=logspace(log10(0.1), log10(1), K);

fprintf('Standard ABC with exchange moves... \n')
rng(321);
params.deltat=fix(1.5*K); params.K = K;
% Start from the posterior (obtained by rejection ABC)
params.theta_in = zeros(2, params.K, C)
for c=1:C
    params.theta_in(:,:,c) = R_rej(randsample(1:length(R_rej), params.K),:)';%repmat(simulations.y', 1, params.K);%
end
fprintf(' every %d local moves for %d seconds \n', params.deltat, params.T)

tic
[Theta_se, X_se, Rej_se, n_se, ne_se, nsw_se, sw_se, TM] = MA_ABC_standard_exchange(C, params, simulations);
fprintf('Exchange moves acceptance rate: %f \n', sw_se'./nsw_se');
b_se = n_se'-ne_se'+1
toc
[~, c_se] = max(b_se(:,1));
fprintf('Standard ABC with exchange moves...done \n\n')

%% Anytime ABC-MCMC algorithm with exchange moves (ABC-APTMC)
% The anytime version of the above algorithm, corrected for length bias.

fprintf('Anytime ABC with exchange moves... \n')
% Ensure median time spent performing local moves before each exchange
% move remains the same as for ABC-PTMC
params.deltat = round(median(TM(TM(:,2,c_se)==1,1,c_se)), 3); params.K = K;
% Start from the posterior (obtained by rejection ABC)
params.theta_in = zeros(2, params.K, C)
for c=1:C
    params.theta_in(:,:,c) = R_rej(randsample(1:length(R_rej), params.K),:)';%repmat(simulations.y', 1, params.K);%
end
fprintf(' every %.2f seconds for %d seconds \n', params.deltat, params.T)

tic
[Theta_a, X_a, Rej_a, n_a, na_e, nsw_a, sw_a, TM_a] = MA_ABC_anytime(C, params, simulations);
fprintf('Exchange moves acceptance rate: %f \n', sw_a'./nsw_a');
b_a = n_a'-na_e'+1 %#ok<*NOPTS>
toc
fprintf('Anytime ABC with exchange moves...done \n')

% Plot the resulting posterior for each chain.
fprintf('Plotting chains... \n')
[~, c_a] = max(b_a(:,1)); 
MA_plot_allchains(Theta_a(:,:,:,c_a), n_a(:, c_a), na_e(:,c_a), theta_true, params, Y, rgb('LightSkyBlue'))
fprintf('Plotting chains...done \n\n')

%% Anytime ABC with exchange moves and no cold updates
fprintf('Anytime ABC with exchange moves and no cold updates... \n')
nocold=1;
tic
[Theta_an, X_an, Rej_an, n_an, ne_an, nsw_an, sw_an, TM_an] = MA_ABC_anytime_all_parallel(C, params, simulations, nocold);
fprintf('Exchange moves acceptance rate: %f \n', sw_an'./nsw_an');
b_an = n_an'-ne_an'+1 %#ok<*NOPTS>
toc

[~, c_an] = max(b_an(:,1));
fprintf('Anytime ABC with exchange moves and no cold updates...done \n\n')
fprintf('Running algorithms...done \n')

fprintf('Saving workspace... ')
save(sprintf('results/MA/MA_comparison_wkspace_%s.mat', datestr(now, 'yyyymmdd_HHMM')))
fprintf('done \n')


%% Performance evaluation
fprintf('Checking performance... \n')
k=1;maxlag=100; Cc = min([sum(b_se(:,1)>maxlag) sum(b_a(:,1)>maxlag) sum(b_an(:,1)>maxlag) sum(b_s(:,1)>maxlag) sum(b_s(:,1)>maxlag)]);
R_s = cell(1, Cc); R_a = cell(1, Cc); R_an = cell(1, Cc); R_se = cell(1, Cc);
s =0; se=0; a=0; an=0;
for c=1:C
    % Vanilla
    if b_s(c)>maxlag
        s=s+1;
        R_s{s} = [Theta_s((ne_s(k, c)+1):n_s(1, c), 1, k, c), Theta_s((ne_s(k, c)+1):n_s(1, c), 2, k, c)];
    end
    % Standard exchange
    if b_se(c,1)>maxlag
        se = se+1;
        R_se{se} = [Theta_se((ne_se(k, c)+1):n_se(1, c), 1, k, c), Theta_se((ne_se(k, c)+1):n_se(1, c), 2, k, c)];
    end
    % Anytime exchange
    if b_a(c,1)>maxlag
        a = a+1;
        R_a{a} = [Theta_a((na_e(k, c)+1):n_a(1, c), 1, k, c), Theta_a((na_e(k, c)+1):n_a(1, c), 2, k, c)];
    end
    % Anytime no cold
    if b_an(c,1)>maxlag
        an = an+1;
        R_an{an} = [Theta_an((ne_an(k, c)+1):n_an(1, c), 1, k, c), Theta_an((ne_an(k, c)+1):n_an(1, c), 2, k, c)];
    end
end

% Efficiency diagnostics
% compute ess and iat
ESS_IAT_e=zeros(1, 4);
ESS_IAT_en=zeros(1, 4);
ESS_IAT_s=zeros(1, 4);
ESS_IAT_se=zeros(1, 4);
maxc=4; 

% Vanilla
[ESS_IAT_s(1), ESS_IAT_s(2)] = autocorr_new(R_s, 1, maxc); [ESS_IAT_s(3), ESS_IAT_s(4)] = autocorr_new(R_s, 2, maxc);
% ABC-PTMC
[ESS_IAT_se(1), ESS_IAT_se(2)] = autocorr_new(R_se, 1, maxc); [ESS_IAT_se(3), ESS_IAT_se(4)] = autocorr_new(R_se, 2, maxc);
% ABC-APTMC
[ESS_IAT_e(1), ESS_IAT_e(2)] = autocorr_new(R_a, 1, maxc); [ESS_IAT_e(3), ESS_IAT_e(4)] = autocorr_new(R_a, 2, maxc);
% ABC-APTMC
[ESS_IAT_en(1), ESS_IAT_en(2)] = autocorr_new(R_an, 1, maxc); [ESS_IAT_en(3), ESS_IAT_en(4)] = autocorr_new(R_an, 2, maxc);

fprintf('Printing performance... \n')
vnames = {'stat', 'ABC_APTMC', 'ABC_APTMC_nocold', 'ABC_PTMC', 'ABC'};
rnames = {'ESS', 'IAT'};
T1 = table(rnames', [ESS_IAT_e(1); ESS_IAT_e(2)],...
	[ESS_IAT_en(1); ESS_IAT_en(2)],...
    	[ESS_IAT_se(1); ESS_IAT_se(2)],...
	[ESS_IAT_s(1); ESS_IAT_s(2)], 'VariableNames', vnames)
T2 = table(rnames', [ESS_IAT_e(3); ESS_IAT_e(4)],...
	[ESS_IAT_en(3); ESS_IAT_en(4)],...
    	[ESS_IAT_se(3); ESS_IAT_se(4)],...
	[ESS_IAT_s(3); ESS_IAT_s(4)], 'VariableNames', vnames)
fprintf('Printing performance...done \n')

% Acf plots
fprintf('Plotting acf... \n')
acf_params.maxlag = maxlag;
acf_params.names = {'Standard ABC', 'ABC-PTMC', 'ABC-APTMC', 'ABC-APTMC-nocold'};
acf_params.linestyle = {'-','-','-','-'};
acf_params.marker = {'*', 'diamond', 'square', 'square'};
acf_params.markersize = {5, 4, 4, 4};
acf_params.col = {rgb('Crimson'), rgb('YellowGreen'), rgb('MediumTurquoise'), rgb('Plum')};
figure;
subplot(2, 1, 1); plot_acf_multi({R_s, R_se, R_a, R_an}, 1, acf_params, 'Sample Autocorrelation function, \theta_1', 1);
subplot(2, 1, 2); plot_acf_multi({R_s, R_se, R_a, R_an}, 2, acf_params, 'Sample Autocorrelation function, \theta_2', 1);
fprintf('Plotting acf...done \n')

fprintf('All done! :D \n')
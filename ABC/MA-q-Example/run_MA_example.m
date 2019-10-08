%% Initialisation
% Applying anytime parallel tempering Monte Carlo (APTMC) to estimate the
% parameters $\theta = \left(\theta_1, \theta_2\right)$ of a MA(q) process.
% Here, $q=2$.

% Various flags to decide which algorithms and diagnostics are run
rejection=0; standard=0; stdexchange = 1;
anytime=1; nocold=1;
acf_plot=1; boxplots=1;
individual_plot=0; three_individual=1;
marginals=0; moremarginals=1;
C=1; %currently only set up for C=1

% Various parameters for simulating observations of the MA(q) process
simulations.q = 2;
theta_true = [0.6; 0.2];
simulations.nsamples = 500;

% Settings of chains for exchange moves
% number of chains for parallel tempering
K=8; params.K = K;
% length bias correction is applied
params.correct=1;

% Burnin and total running time of algorithms
b=100; params.T = b+300; params.burnin = b;
params.N=1e7;


%%
% Define the radius $\varepsilon$ of the ball associated with each chain. It can be defined so that each chain
% starts with a larger epsilon which pogressively shrinks during the burnin
% period.
epsilone = logspace(log10(0.02), log10(1), K); %logspace(log10(2), log10(2.98), K);
epsilons = logspace(log10(0.02), log10(1), K);
params.epsilon = zeros(params.T, K);
for k=1:K
    params.epsilon(1:b,k) = linspace(epsilone(k), epsilons(k), b);
    params.epsilon((b+1):params.T,k) = epsilons(k);
end

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
fprintf('Simulating observations...done \n')

%% ALGORITHMS runs
%% ABC Rejection sampling
% For reference, obtain an approximation of the posterior by generating
% $N_r$ samples from the prior, simulating observations and only retaining
% those which hit the ball of radius $\varepsilon = 0.02$.
rng('default');
if(rejection)          
    fprintf('Rejection ABC... \n') %#ok<*UNRCH>
    N_r = 1e7;
    % run on a parallel pool to speed up algorithm
    if(isempty(gcp('nocreate')))
        poolobj = parpool(min(N_r, 4));
    end
    tic
    % Rejection ABC
    [R_rej, X_r, Rej_r, n_r] = MA_ABC_rejection(N_r, simulations, epsilons(1));
    toc
    dlmwrite(sprintf('results/ABC/MA/MA_ABC_rejection_%d.csv', n_r), R_rej);
    fprintf('Acceptance rate = %f \n\n', n_r/N_r)
    fprintf('Rejection ABC...done \n\n')
    delete(poolobj)
end
R_rej = dlmread(sprintf('results/ABC/MA/MA_ABC_rejection_%d.csv', 13189));
fprintf('Rejection ABC: Acceptance rate = %f percent \n', 13189/1e7*100)

%%
% This algorithm has an extremely low acceptance rate and is therefore
% highly inefficient, which is why it is a good idea to resort to MCMC.

%% Standard ABC-MCMC algorithm
% The vanilla MCMC algorithm which employs ABC in its local moves. The
% algorithm runs for 3600 seconds on a chain associated with the ball of
% radius $\varepsilon = 0.02$.

if(standard)
    fprintf('Standard ABC... \n') %#ok<*UNRCH>
    params.deltat = params.T+1000; params.K = 1;
    % Start from the posterior (obtained by rejection ABC)
    params.theta_in = R_rej(randsample(1:length(R_rej), params.K),:)';
    
    % Standard deviation for Gaussian random walk proposal for local moves
    params.rho=0.1;
    
    tic
    [Theta, X, Rej, n, ne] = MA_ABC_standard(C, params, simulations, params.epsilon(:,1));
    b_s = n'-ne'+1
    toc
    fprintf('Standard ABC...done \n\n')
    % Select run which returned the most samples
    [~, c] = max(n(1,:)-ne(1,:)+1);
    R_s = Theta(ne(1,c):n(1,c), :, 1, c);
    
    % Break down the composition of the resulting chain and compute local
    % moves acceptance rate (as percentage)
    fprintf(' \n --- CHAIN 1, epsilon = %.2f --- \n \n', params.epsilon(1))
    fprintf('Chain is made up of: \n')
    fprintf('Rejection due to prior: %f percent |',sum(Rej(ne:n)==2)/b_s*100)
    fprintf(' due to race: %f percent \n',sum(Rej(ne:n)==1)/b_s*100)
    fprintf('Accepted local move: %f percent \n',sum(Rej(ne:n)==0)/b_s*100)
    
    % Record cold chain
    dlmwrite(sprintf('results/ABC/MA/MA_ABC_standard_%d.csv', b_s), R_s);
    % R_s = dlmread(sprintf('results/ABC/MA/MA_ABC_standard_%d.csv', 33634));
    
end

%% Standard ABC-MCMC algorithm with exchange moves (ABC-PTMC)
% The above algorithm but run on $K=8$ chains, each targetting a posterior
% associated with a different $\varepsilon^k$. The first, _cold_ chain is the one
% associated with ball of radius $\varepsilon^1 = 0.02$, the _warmer_ chains
% target posteriors associated with balls of larger radii. Every $\delta_t$
% local moves an exchange move occurs between two uniformly randomly chosen adjacent chains.

rng(678451)
if(stdexchange)
    
    % Standard deviation for Gaussian random walk proposal for local moves
    params.rho=logspace(log10(0.1), log10(1), K);
    
    fprintf('Standard ABC with exchange moves... \n')
    params.deltat=fix(K/2); params.K = K;
    % Start from the posterior (obtained by rejection ABC)
    params.theta_in = R_rej(randsample(1:length(R_rej), params.K),:)';
    fprintf(' every %d local moves for %d seconds \n', params.deltat, params.T)
    
    tic
    [Theta_se, X_se, Rej_se, n_se, ne_se, nsw_se, sw_se, TM] = MA_ABC_standard_exchange(C, params, simulations);
    fprintf('Exchange moves acceptance rate: %f \n', sw_se'./nsw_se');
    b_se = n_se'-ne_se'+1 %#ok<*NOPTS>
    toc
    fprintf('Standard ABC with exchange moves...done \n\n')
    
    % Break down the composition of each chain and compute local and
    % exchange moves acceptance rates (as percentages)
    chains_s.local_rejected_prior = zeros(1,K); chains_s.local_rejected_race = zeros(1,K);
    chains_s.local_accepted = zeros(1,K);
    chains_s.exchange_rejected = zeros(1,K); chains_s.exchange_accepted = zeros(1,K);
    rates_s.local = zeros(1,K); rates_s.exchange = zeros(1,K);
    for k=1:K
        chains_s.local_rejected_prior(k) = sum(Rej_se(ne_se(k):n_se(k), k)==2)/b_se(k)*100;
        chains_s.local_rejected_race(k) = sum(Rej_se(ne_se(k):n_se(k), k)==1)/b_se(k)*100;
        chains_s.local_accepted(k) = sum(Rej_se(ne_se(k):n_se(k), k)==0)/b_se(k)*100;
        chains_s.exchange_rejected(k) = sum(Rej_se(ne_se(k):n_se(k), k)==3)/b_se(k)*100;
        chains_s.exchange_accepted(k) = sum(Rej_se(ne_se(k):n_se(k), k)==4)/b_se(k)*100;
         rates_s.local(k) = sum(Rej_se(ne_se(k):n_se(k), k)==0)/sum(Rej_se(ne_se(k):n_se(k), k)<3)*100 ;
        rates_s.exchange(k) = sum(Rej_se(ne_se(k):n_se(k), k)==4)/sum(Rej_se(ne_se(k):n_se(k), k)>=3)*100 ;
    end
    chains_s
    rates_s
    
    % Plot the resulting posterior for each chain.
    fprintf('Plotting chains... \n')
    [~, c] = max(n_se(1,:)-ne_se(1,:)+1);
    MA_plot_allchains(Theta_se(:,:,:,c), n_se(:,c), ne_se(:,c), theta_true, params, Y, rgb('YellowGreen'))
    fprintf('Plotting chains...done \n\n')
    
    % Record cold chain
    R_se = Theta_se(ne_se(1,c):n_se(1,c), :, 1, c);
    dlmwrite(sprintf('results/ABC/MA/MA_ABC_stexchange_%d.csv', b_se(1)), R_se);
    % R_se = dlmread(sprintf('results/ABC/MA/MA_ABC_stdexchange_%d.csv', 163491));    
    
    
    % Record timeline of local and exchange moves
    [~, i] = max(cumsum(TM(:,1))>b); TM_t = TM(i:end,:);
    dlmwrite(sprintf('results/ABC/MA/MA_ABC_stdexchange_TM_%d.csv', b_se(1)), TM_t);
    % TM_t = dlmread(sprintf('results/ABC/MA/MA_ABC_stdexchange_TM_%d.csv', 163491));
    
end

%% Anytime ABC-MCMC algorithm with exchange moves (ABC-APTMC)
% The anytime version of the above algorithm, corrected for length bias.
% Again run on $K=8$ chains, each targetting a posterior associated with a
% different $\varepsilon^k$. This time, exchange moves occur every
% $\delta^a_t$ seconds, where now $\delta^a_t$ is the median time the
% previous algorithm took to perform $\delta_t$ local moves.

rng(23672);
if(anytime)
    fprintf('Anytime ABC with exchange moves... \n')
    % Ensure median time spent performing local moves before each exchange
    % move remains the same as for ABC-PTMC
    params.deltat = round(median(TM_t(TM_t(:,2)==1,1)), 3); params.K = K;
    % Start from the posterior (obtained by rejection ABC)
    params.theta_in = R_rej(randsample(1:length(R_rej), params.K),:)';
    fprintf(' every %.2f seconds for %d seconds \n', params.deltat, params.T)
    
    tic
    [Theta_e, X_e, Rej_e, n_e, ne_e, nsw_e, sw_e, TM_a] = MA_ABC_anytime_all(C, params, simulations);
    fprintf('Exchange moves acceptance rate: %f \n', sw_e'./nsw_e');
    b_e = n_e'-ne_e'+1 %#ok<*NOPTS>
    toc
    fprintf('Anytime ABC with exchange moves...done \n\n')
    
    % Break down the composition of each chain and compute local and
    % exchange moves acceptance rates (as percentages)
    chains_a.local_rejected_prior = zeros(1,K); chains_a.local_rejected_race = zeros(1,K);
    chains_a.local_accepted = zeros(1,K);
    chains_a.exchange_rejected = zeros(1,K); chains_a.exchange_accepted = zeros(1,K);
    rates_a.local = zeros(1,K); rates_a.exchange = zeros(1,K);
    for k=1:K        
        chains_a.local_rejected_prior(k) = sum(Rej_e(ne_e(k):n_e(k), k)==2)/b_e(k)*100;
        chains_a.local_rejected_race(k) = sum(Rej_e(ne_e(k):n_e(k), k)==1)/b_e(k)*100;
        chains_a.local_accepted(k) = sum(Rej_e(ne_e(k):n_e(k), k)==0)/b_e(k)*100;
        chains_a.exchange_rejected(k) = sum(Rej_e(ne_e(k):n_e(k), k)==3)/b_e(k)*100;
        chains_a.exchange_accepted(k) = sum(Rej_e(ne_e(k):n_e(k), k)==4)/b_e(k)*100;
        rates_a.local(k) = sum(Rej_e(ne_e(k):n_e(k), k)==0)/sum(Rej_e(ne_e(k):n_e(k), k)<3)*100 ;
        rates_a.exchange(k) = sum(Rej_e(ne_e(k):n_e(k), k)==4)/sum(Rej_e(ne_e(k):n_e(k), k)>=3)*100;        
    end
    chains_a
    rates_a
    
    % Plot the resulting posterior for each chain.
    fprintf('Plotting chains... \n')
    [~, c] = max(n_e(1,:)-ne_e(1,:)+1);
    MA_plot_allchains(Theta_e(:,:,:,c), n_e(:,c), ne_e(:,c), theta_true, params, Y, rgb('LightSkyBlue'))
    fprintf('Plotting chains...done \n\n')
    
    % Record cold chain
    R_e = Theta_e(ne_e(1,c):n_e(1,c), :, 1, c);
    dlmwrite(sprintf('results/ABC/MA/MA_ABC_anytime_%d.csv', b_e(1)), R_e);
    % R_e = dlmread(sprintf('results/ABC/MA/MA_ABC_anytime_%d.csv', 149316));    
    
    % Record timeline of local and exchange moves
    [~, j] = max(cumsum(TM_a(:,1))>b); TM_a_t = TM_a(j:end,:);
    dlmwrite(sprintf('results/ABC/MA/MA_ABC_anytime_TM_%d.csv', b_e(1)), TM_a_t);
    % TM_a_t = dlmread(sprintf('results/ABC/MA/MA_ABC_anytime_TM_%d.csv', 149316));
   
end


%% Anytime ABC-MCMC algorithm with exchange moves and no cold local moves (ABC-APTMC-nocold)
% The same algorithm as above, except local moves do not occur on the cold
% chain, which is therefore entirely made up of samples from exchange moves.

rng(396877);
if(nocold)
    fprintf('Anytime ABC with exchange moves and no cold updates... \n')
    params.deltat = round(median(TM_t(TM_t(:,2)==1,1)), 3); params.K = K;
    % Start from the posterior (obtained by rejection ABC)
    params.theta_in = R_rej(randsample(1:length(R_rej), params.K),:)';
    fprintf(' every %.2f seconds for %d seconds \n', params.deltat, params.T)
    % No local moves on the cold chain.
    nocold = 1;
    
    tic
    [Theta_en, X_en, Rej_en, n_en, ne_en, nsw_en, sw_en, TM_an] = MA_ABC_anytime_all(C, params, simulations, nocold);
    fprintf('Exchange moves acceptance rate: %f \n', sw_en'./nsw_en');
    b_en = n_en'-ne_en'+1 %#ok<*NOPTS>
    toc
    fprintf('Anytime ABC with exchange moves and no cold updates...done \n\n')
    
    % Break down the composition of each chain and compute local and
    % exchange moves acceptance rates (as percentages)
    chains_an.local_rejected_prior = zeros(1,K); chains_an.local_rejected_race = zeros(1,K);
    chains_an.local_accepted = zeros(1,K);
    chains_an.exchange_rejected = zeros(1,K); chains_an.exchange_accepted = zeros(1,K);
    rates_an.local = zeros(1,K); rates_an.exchange = zeros(1,K);
    for k=1:K
        chains_an.local_rejected_prior(k) = sum(Rej_en(ne_en(k):n_en(k), k)==2)/b_en(k)*100;
        chains_an.local_rejected_race(k) = sum(Rej_en(ne_en(k):n_en(k), k)==1)/b_en(k)*100;
        chains_an.local_accepted(k) = sum(Rej_en(ne_en(k):n_en(k), k)==0)/b_en(k)*100;
        chains_an.exchange_rejected(k) = sum(Rej_en(ne_en(k):n_en(k), k)==3)/b_en(k)*100;
        chains_an.exchange_accepted(k) = sum(Rej_en(ne_en(k):n_en(k), k)==4)/b_en(k)*100;
        rates_an.local(k) = sum(Rej_en(ne_en(k):n_en(k), k)==0)/sum(Rej_en(ne_en(k):n_en(k), k)<3)*100 ;
        rates_an.exchange(k) = sum(Rej_en(ne_en(k):n_en(k), k)==4)/sum(Rej_en(ne_en(k):n_en(k), k)>=3)*100 ;
    end
    chains_an
    rates_an
    
    % Plot the resulting posterior for each chain.
    fprintf('Plotting chains... \n')
    [~, c] = max(n_en(1,:)-ne_en(1,:)+1);
    MA_plot_allchains(Theta_en(:,:,:,c), n_en(:,c), ne_en(:,c), theta_true, params, Y, rgb('Plum'))
    fprintf('Plotting chains...done \n\n')
    
    % Record cold chain
    R_en = Theta_en(ne_en(1,c):n_en(1,c), :, 1, c);
    dlmwrite(sprintf('results/ABC/MA/MA_ABC_anytime_nocold_%d.csv', b_en(1)), R_en);
    % R_en = dlmread(sprintf('results/ABC/MA/MA_ABC_anytime_nocold_%d.csv', 171205));    

    % Record timeline of local and exchange moves
    [~, j] = max(cumsum(TM_an(:,1))>b); TM_an_t = TM_an(j:end,:);
end

%% DIAGNOSTICS
% In this section, diagnostics are performed on the outputs of all ABC-MCMC algorithms.

%% Performance comparison
% All ABC-MCMC algorithms ran for the same amount of time. To compare their
% efficiency, sample autocorelation function plots are drawn and the
% effective sample size (ESS) and integrated autocorrelation time (iat) are
% computed using the  initial monotone positive sequence estimator (IMSE)

if(acf_plot)
    % Plot the sample acf for the ABC, ABC-PTMC and ABT-APTMC algorithms
    maxlag=min([201, b_se(1), b_e(1), b_s, b_en(1)])-1;
    figure;
    subplot(2, 1, 1); plot_acf_four(R_s(:,1), R_se(:,1), R_e(:,1), R_en(:,1),...
        maxlag, 'Sample Autocorrelation Function, \theta_1'); 
    subplot(2, 1, 2); plot_acf_four(R_s(:,2), R_se(:,2), R_e(:,2), R_en(:,2),...
        maxlag, 'Sample Autocorrelation Function, \theta_2'); 
    
    % Compute ess and iat for all ABC-MCMC algorithms
    imse=1;
    % Standard ABC
    [ess.standard1, iat.standard1] = ESS_IAT(R_s(:,1), maxlag, imse);
    [ess.standard2, iat.standard2] = ESS_IAT(R_s(:,2), maxlag, imse);
    % ABC-PTMC
    [ess.stdexchange1, iat.stdexchange1] = ESS_IAT(R_se(:,1), maxlag, imse);
    [ess.stdexchange2, iat.stdexchange2] = ESS_IAT(R_se(:,2), maxlag, imse);
    % ABC-APTMC
    [ess.anytime1, iat.anytime1] = ESS_IAT(R_e(:,1), maxlag, imse);
    [ess.anytime2, iat.anytime2] = ESS_IAT(R_e(:,2), maxlag, imse);
    % ABC-APTMC-nocold
    [ess.anytimen1, iat.anytimen1] = ESS_IAT(R_en(:,1), maxlag, imse);
    [ess.anytimen2, iat.anytimen2] = ESS_IAT(R_en(:,2), maxlag, imse)
end

%% Plots of marginals
% For reference, the marginal posteriors for $\theta_1$ and $\theta_2$ associated with the
% ball of radius $\varepsilon=0.02$ can approximated using a kernel density estimate of the ABC
% rejection sampling output. Kernel density estimates are also be obtained
% for all the other algorithms.
if(moremarginals)
    fprintf('Plotting marginals... \n')
    
    % True marginals can also be computed in this example, but it takes some time.
%     x1 = linspace(-2, 2, 1000);
%     y2 = linspace(-1, 1, 1000);
%     
%     %F = @(x) MA_likelihoodx(x, y2, Y, 1);
%     Zx = 1.23084894785056e-113; %normalisation constant
%     %Zx = integral(F, -2, 2);
%     %     F = @(x) MA_likelihoodx(x, y2, Y, Zx);
%     %     integral(F, -2, 2)
%     Stx = MA_likelihoodx(x1, y2, Y, Zx);
%     
%     
%     %G = @(y) MA_likelihoody(x1, y, Y, 1);
%     Zy = 4.369447225175966e-23; %normalisation constant
%     %Zy = integral(G, -1, 1);
%     %     G = @(y) MA_likelihoody(x1, y, Y, Zy);
%     %     integral(G, -1, 1)
%     Sty = MA_likelihoody(x1, y2, Y, Zy);
    
    figure;
    %THETA_1 kernel density estimates
    [~, xdensity_r, xmesh_r] = kde(R_rej(:,1), 2^10);
    S_s= cell2mat(R_s'); [~, xdensity, xmesh] = kde(S_s(:,1), 2^10);
    S_e= cell2mat(R_e'); [~, xdensity_e, xmesh_e] = kde(S_e(:,1), 2^10);
    S_se= cell2mat(R_se'); [~, xdensity_se, xmesh_se] = kde(S_se(:,1), 2^10);
    S_en= cell2mat(R_en'); [~, xdensity_en, xmesh_en] = kde(S_en(:,1), 2^10); 
       
    
    subplot(1, 2, 1); hold on;
    % True value
    plot([0.6 0.6], [0 3.5], ...
        'linestyle', '--', ...
        'color', rgb('Violet'));
    % True marginal posterior
    plot(x1, Stx, 'Color', 'black', 'LineStyle', '-.');
    % Rejection ABC posterior
    plot(xmesh_r, xdensity_r, 'Color', 'black', 'LineStyle', '-.', 'LineWidth', 2); xlim([0 1]);
    % Standard ABC posterior
    plot(xmesh, xdensity, 'Color', rgb('Crimson'));
    % ABC-APTMC posterior
    plot(xmesh_se, xdensity_se, 'Color', rgb('YellowGreen'), 'LineWidth', 1);
    % ABC-APTMC posterior
    plot(xmesh_e, xdensity_e, 'Color', rgb('MediumTurquoise'), 'LineWidth', 1);
    % ABC-APTMC-nocold posterior
    plot(xmesh_en, xdensity_en, 'Color', rgb('Plum'), 'LineWidth', .5);  xlabel('\theta_1');    
    
    %THETA_2 kernel density estimates
    [~, ydensity_r, ymesh_r] = kde(R_rej(:,2), 2^10);
    [~, ydensity, ymesh] = kde(S_s(:,2), 2^10);
    [~, ydensity_e, ymesh_e] = kde(S_e(:,2), 2^10);
    [~, ydensity_se, ymesh_se] = kde(S_se(:,2), 2^10);
    [~, ydensity_en, ymesh_en] = kde(S_en(:,2), 2^10);
    clear('S_s', 'S_e', 'S_se', 'S_en')
    subplot(1, 2, 2); hold on;
    % True value
    plot([0.2 0.2], [0 2], ...
        'linestyle', '--', ...
        'color', rgb('Violet'));
    % True marginal posterior
    plot(y2, Sty, 'Color', 'black', 'LineStyle', '-.');
    % Rejection ABC posterior
    plot(ymesh_r, ydensity_r, 'Color', 'black', 'LineStyle', '-.', 'LineWidth', 2);xlim([-1 1.1]);
    % Standard ABC posterior
    plot(ymesh, ydensity, 'Color', rgb('Crimson'));
    % ABC-PTMC posterior
    plot(ymesh_se, ydensity_se, 'Color', rgb('YellowGreen'), 'LineWidth', 1);
    % ABC-APTMC posterior
    plot(ymesh_e, ydensity_e, 'Color', rgb('MediumTurquoise'), 'LineWidth', 1);
    % ABC-APTMC no cold updates posterior
    plot(ymesh_en, ydensity_en, 'Color', rgb('Plum'), 'LineWidth', .5); xlabel('\theta_2');       
    
    fprintf('Plotting marginals...done \n\n')
end

%% Plot a single chain
% Used here for displaying the posterior returned by single chain in the
% standard ABC algorithm.

if(individual_plot)
    k=1;
    S1 = Theta((ne(k, c)+1):n(1, c), 1, k, c);
    S2 = Theta((ne(k, c)+1):n(1, c), 2, k, c);
    p = gkde2([S1,S2]);
    
    % True posterior can be computed but takes some time
    %     [cx, cy] = fixxy(p.x, p.y);
    %
    %     F = @(t1, t2) MA_likelihood2(t1, t2, Y);
    %     Z = integral2(F, -2, 2, -1, 1);
    %     F = @(t1, t2) Z^(-1)*MA_likelihood2(t1, t2, Y);
    %     Z = Z*integral2(F, -2, 2, -1, 1);
    %     S = MA_likelihood2(cx, cy, Y, Z);
    
    % Prior 'triangle'
    xx = [-2 0 2]; yy = [1 -1 1];
    figure;  fill(xx,yy,rgb('Lavender'))
    % Scatter plot of individual posterior
    hold on; scatter(S1, S2, 10, rgb('Crimson'), 'filled');
    % True posterior    
    %     zlevs = get_levels(S, 5)
    %     contour(cx, cy, S, 'LineColor', 'black', 'LineStyle', '-.', 'LevelList', zlevs, 'LineWidth', .1);
    % Contour plot of 2D kde estimate of the individual posterior
    p.zlevs = get_levels(p.pdf, 8);
    contour(p.x,p.y,p.pdf, 'LevelList', p.zlevs);
    % True value
    scatter(theta_true(1), theta_true(2), 50, rgb('LawnGreen'), 'filled', 'pentagram');
    xlabel('\theta_1'); ylabel('\theta_2');    
end

%% Display multiple cold chains on a single plot
% Used here to compare the cold posterior for the ABC-MCMC algorithms.

if(three_individual)        
    % 2D kde estimate of the ABC posterior using the output from the
    % rejection ABC algorithm
    p = gkde2(R_rej);
    p.zlevs = get_levels(p.pdf, 8);
    
    figure;  
    % Prior 'triangle'    
    xx = [-2 0 2]; yy = [1 -1 1];
    fill(xx,yy,rgb('Lavender'))
    hold on;
    % Standard ABC
    h1=scatter(R_s{c_s}(:,1), R_s{c_s}(:,2), 10, rgb('Crimson'), 'filled',...
        'DisplayName', 'Standard ABC');
    % ABC-PTMC
    h2=scatter(R_se{c_se}(:,1), R_se{c_se}(:,2), 7, rgb('YellowGreen'), 'filled',...
        'DisplayName', 'ABC-PTMC');
    % ABC-APTMC
    h3=scatter(R_e{c_e}(:,1), R_e{c_e}(:,2), 5, rgb('MediumTurquoise'), 'filled',...
        'DisplayName', 'ABC-APTMC');
    %ABC-APTMC-nocold
    h4=scatter(R_en{c_en}(:,1), R_en{c_en}(:,2), 3, rgb('Plum'), 'filled',...
        'DisplayName', 'ABC-APTMC-nocold');
    % Contour plot of 2D kde estimate of the ABC posterior
    contour(p.x, p.y, p.pdf, 'LevelList', p.zlevs);
    % True value
    scatter(theta_true(1), theta_true(2), 75, rgb('Gold'), 'filled',...
        'pentagram', 'MarkerEdgeColor',rgb('OrangeRed'));
    
    lgd = legend([h1 h2 h3 h4]);
    lgd.Location='southwest';
    xlabel('\theta_1');
    ylabel('\theta_2');
end

%% Time plots
% Two kinds of plots are returned here:
%
% # Timelines of local and exchange moves for the ABC-PTMC and ABC-APTMC
% algorithms for the first 0.1 seconds of their runs, to compare the time
% each algorithm spends performing each type of update.
% # Boxplots of all the times taken to complete local moves (since
% displaying them all on a timeline would make it unreadable).
% 

if(boxplots)     
    % Timeline plots
    [~, i] = max(cumsum(TM(:,1))>b); TM_t = TM(i:end,:);
    [~, j] = max(cumsum(TM_a(:,1))>b); TM_a_t = TM_a(j:end,:);
    timeplots(TM_t, 'Standard', rgb('MediumBlue'), [0, 1])
    timeplots(TM_a_t, 'Anytime', rgb('DarkOrange'), [0, 1])
    
    % Boxplots of local move times
    % Create ‘NaN’ Array
    C3 = NaN(max([sum(TM_t(:,2)==1), sum(TM_a_t(:,2)==1)]), 2); 
    C3(1:sum(TM_t(:,2)==1),1) = TM_t(TM_t(:,2)==1,1);
    C3(1:sum(TM_a_t(:,2)==1), 2) =  TM_a_t(TM_a_t(:,2)==1,1);
    figure; bxplt=boxplot(C3, {'ABC-PTMC', 'ABC-APTMC'});
    ylabel('local move times (seconds)');
    bxplt.YLim = [0.0 0.1];
    title('Boxplot of time taken to perform local moves');
end

fprintf('Saving workspace... \n')
save(sprintf('results/MA/MA_example_wkspace_%s.mat', datestr(now, 'yyyymmdd_HHMM')))
fprintf('Saving workspace...done \n')


fprintf('All done! :D \n')
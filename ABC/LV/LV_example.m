%% Single Processor Example
% Applying approximate Bayesian computation (ABC) with parallel tempering
% Monte Carlo (PTMC) and anytime parallel tempering Monte Carlo (APTMC) to
% estimate the parameters $\theta = \left(\theta_1, \theta_2,
% \theta_3\right)$ of a stochastic Lotka-Volterra predator-prey model.

% Various flags to decide which algorithms and diagnostics are run
 clear;
run_rejection = 0; run_vanilla = 1; run_stdexchange = 1; run_anytime=1;
run_one = 1; run_multi = 0;
% Plots and diagnostics
acf_plots = 1; time_plots = 0; density_plots = 0; 
ess_plot = 1; chain_plots = 0; marginal_plots = 0;

LV_initialisation
fprintf('I will be done on:')
disp(datetime(clock + [0 0 0 0 0 (run_vanilla+run_anytime+run_stdexchange)*params.T]))

%% Standard ABC-MCMC algorithm
% The vanilla MCMC algorithm which employs ABC in its local moves. The
% algorithm runs on a single chain associated with the ball of
% radius $\varepsilon = 1$.

if(run_vanilla)
    % Run LV_standard only performing vanilla ABC-MCMC
    rejection = 0; vanilla = 1; 
    LV_standard    
end

%% Standard ABC-MCMC algorithm with exchange moves (ABC-PTMC-1)
% The above algorithm but run sequentially on $K=6$ chains, each targetting a posterior
% associated with a different $\varepsilon^k$. The first, _cold_ chain is the
% one associated with ball of radius $\varepsilon^1 = 1$, the _warmer_
% chains target posteriors associated with balls of larger radii. Every
% $\delta_t$ local moves an exchange move occurs between two uniformly
% randomly chosen adjacent chains.

if(run_stdexchange)
    % Run LV_stdexchange only performing single processor ABC-PTMC-1
    run_single_stdexchange = 1; run_multi_stdexchange = 0;
    LV_stdexchange
end

%% Anytime ABC-MCMC algorithm with exchange moves (ABC-APTMC-1)
% The anytime version of the ABC_APTMC-1 algorithm, corrected for length bias.
% Again run on $K=6$ chains, each targetting a posterior associated with a
% different $\varepsilon^k$. This time, exchange moves occur every
% $\delta^a_t$ seconds, where now $\delta^a_t$ is the median time the
% previous algorithm took to perform $\delta_t$ local moves.

if(run_anytime)
    % Run LV_anytime only performing single processor ABC-APTMC-1
    run_single_anytime = 1; run_multi_anytime = 0;
    LV_anytime
end

%% Performance comparison
% All ABC-MCMC algorithms ran for the same amount of time. To compare their
% efficiency, sample autocorelation function plots are drawn and the
% effective sample size (ESS) and integrated autocorrelation time (iat) are
% computed using the initial monotone positive sequence estimator (IMSE).

if(acf_plots)
    % Compute ess and iat for all ABC-MCMC algorithms    
    maxc = 6;    
    
    % Standard ABC
    [ess.standard, iat.standard] = autocorr_new_multi(S, maxc);
    % ABC-PTMC-1
    [ess.stdexchange, iat.stdexchange] = autocorr_new_multi(S_e, maxc);    
    % ABC-APTMC-1
    [ess.anytime, iat.anytime] = autocorr_new_multi(A_e, maxc);            
               
    fprintf("ESS:\n"); disp(ess);
    fprintf("IAT:\n"); disp(iat);       
    
    fprintf('Standard ABC has IAT: \n') ;
    fprintf(' %.3f larger on average than ABC-PTMC-1 \n', mean(iat.standard./iat.stdexchange));
    fprintf(' %.3f larger on average than ABC-APTMC-1 \n', mean(iat.standard./iat.anytime));
    
    fprintf('Standard ABC has ESS: \n') ;
    fprintf(' %.3f smaller on average than ABC-PTMC-1 \n', mean(ess.stdexchange./ess.standard));
    fprintf(' %.3f smaller on average than ABC-APTMC-1 \n', mean(ess.anytime./ess.standard));   
    
    % ACF plots
    acf_params.maxlag = 200;
    acf_params.names = {'Standard ABC', 'ABC-PTMC', 'ABC-APTMC'};
    acf_params.linestyle = {':','--',':'};
    acf_params.marker = {'diamond', 'square', 'diamond'};
    acf_params.markersize = {4, 3, 2, 2};
    acf_params.col = { [12 195 82] ./ 255, rgb('DarkBlue'),  rgb('DarkOrange')};
    [R_out] = acf_keep_maxlag({S, S_e, A_e}, {b_m(1,:), b_e(1,:),b_a(1, :)}, acf_params.maxlag);
    figure;
    for i=1:3
        subplot(3, 1, i); plot_acf_multi(R_out, i, acf_params, sprintf('\\theta_%d', i), 1);
    end   
end

%% Chain plots
% Plotting the chains output by the various ABC-MCMC algorithms
if(chain_plots)
    % Single processor chains
    plot_multi_chains(S, true_theta, 'Standard ABC on K processors') %#ok<*UNRCH>
    plot_multi_chains(S_e, true_theta, 'ABC-PTMC-1 on K processors')
    plot_multi_chains(A_e, true_theta, 'ABC-APTMC-1 on K processors')    
end

%% Marginal plots
% Plot histograms of the marginal posteriors output from ABC-MCMC
% algorithms and compare them to the ABC rejection sampling reference
if(marginal_plots)
    nbins = 70;
    plot_marginal_posteriors(cell2mat(S'), S_rej, nbins, 'Standard ABC on K processors')
    plot_marginal_posteriors(cell2mat(S_e'), S_rej, nbins, 'ABC-PTMC-1 on K processors')
    plot_marginal_posteriors(cell2mat(A_e'), S_rej, nbins, 'ABC-APTMC-1 on K processors')    
end

if(density_plots)
    % Display all marginal densities on the same plot
    den = [2^8 2^16 2^8
        2^8 2^8 2^8
        2^8 2^8 2^8
        2^8 2^8 2^8];
        plot_multi_densities({cell2mat(S_e'), cell2mat(A_e'), cell2mat(S')},...
            S_rej, den);    
end

%% Timeplots
% Timelines of local and exchange moves for the ABC-(A)PTMC algorithms for the
% first 300 seconds of their runs, to observe the time the algorithms spend
% performing each type of update.

if(time_plots)
    % Timeline plots (one processor)
    time_interval = [0 min(300, params.T)]; %[max(0, params.T-300) params.T];
    w = 1; TMS = TM{w};
    TMA = TM_a{w};
    
    figure;
    subplot(2, 1, 1);
    timeplots(TMS, 'ABC-PTMC-1: Timeline of local and exchange moves',...
        {rgb('MediumBlue'), rgb('LightCyan')}, {rgb('LightCyan'), rgb('MediumBlue')}, time_interval);
    subplot(2, 1, 2);
    timeplots(TMA, 'ABC-APTMC-1: Timeline of local and exchange moves',...
        {rgb('DarkOrange'), rgb('SandyBrown')}, {rgb('DarkRed'), rgb('DarkOrange')}, time_interval);
end

%% ESS over time
% Plot of effective sample size of all algorithms over time. Based on
% sample sizes taken at regular time intervals and final iats.
% A linear regression is fitted to each algorithm for visualisation.
% Confidence intervals can also be drawn. 
if ess_plot   
    plot_essovertime(W, {TM_a, TM, TM_m}, iat, T_s, 0);
end

fprintf('Saving workspace... ')
save(sprintf('results/ABC/LV/LV_example_wkspace_%s.mat', datestr(now, 'yyyymmdd_HHMM')))
fprintf('done \n')

delete(gcp);
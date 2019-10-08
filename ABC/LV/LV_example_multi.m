%% Initialisation
% Applying approximate Bayesian computation (ABC) with parallel tempering
% Monte Carlo (PTMC) and anytime parallel tempering Monte Carlo (APTMC) to
% estimate the parameters $\theta = \left(\theta_1, \theta_2,
% \theta_3\right)$ of a stochastic Lotka-Volterra predator-prey model.

% Various flags to decide which algorithms and diagnostics are run
run_rejection = 0; run_vanilla =0; run_stdexchange = 1; run_anytime=1;
run_one = 0; single_on_multi = 0;
run_multi = 1; multi_twice = 1;
% Plots and diagnostics
acf_plots = 0; time_plots = 0; box_plots = 0;
density_plots = 0; mean_plots = 0;
chain_plots = 0; marginal_plots = 0; all_chains_ess=0;


LV_initialisation_multi
fprintf('I will be done on:')
disp(datetime(clock + [0 0 0 0 0 (2*run_one*ceil(ntimes/4)+(run_anytime+run_stdexchange)*ntimes)*params.T]))

%% Standard ABC-MCMC algorithm
% The vanilla MCMC algorithm which employs ABC in its local moves. The
% algorithm runs for on a chain associated with the ball of radius $\varepsilon = 1$.

if(run_vanilla)
    % Run LV_standard only performing Single processor ABC-MCMC
    standard_rejection=0; 
    run_single_standard=1-single_on_multi; 
    run_multi_standard=single_on_multi;
    LV_standard    
end
     
    
%% Multi-processor standard ABC-MCMC algorithm with exchange moves (ABC-PTMC-W)
% The above algorithm but run in parallel on 4 workers, with $Kk=2$
% chains per processor (i.e. K=8 chains in total), each chain targets a
% posterior associated with a different $\varepsilon^k$ and exchange moves
% are performed between workers.

if(run_stdexchange)    
    % Run LV_stdexchange only performing multi-processor ABC-PTMC-K
    run_single_stdexchange=0; run_single_stdexchange_on_multi=0;
    run_multi_stdexchange=1-multi_twice; run_multi_stdexchange_multi=multi_twice;
    LV_stdexchange
end

%% Multi-processor anytime ABC-MCMC algorithm with exchange moves (ABC-APTMC-W)
% The above algorithm but run in parallel on 4 workers, with $Kk=2$
% chains per processor (i.e. K=8 chains in total), each chain targets a
% posterior associated with a different $\varepsilon^k$ and exchange moves
% are performed between workers.

if(run_anytime)
    % Run LV_anytime only performing multi-processor ABC-APTMC-K
    
    run_single_anytime=0; run_single_anytime_on_multi=0;
    run_multi_anytime=1-multi_twice; run_multi_anytime_multi=multi_twice;
    LV_anytime
    
end

if run_one
    W = ntimes;
    % Run LV_stdexchange only performing single-processor ABC-PTMC-1
    run_single_stdexchange=0; run_single_stdexchange_on_multi=1;
    run_multi_stdexchange=0; run_multi_stdexchange_multi=0;
    LV_stdexchange
    
    % Run LV_anytime only performing multi-processor ABC-APTMC-1    
    run_single_anytime=0; run_single_anytime_on_multi=1;
    run_multi_anytime=0; run_multi_anytime_multi=0;
    LV_anytime
    
end

%% Performance comparison
% All ABC-MCMC algorithms ran for the same amount of time. To compare their
% efficiency, sample autocorelation function plots are drawn and the
% effective sample size (ESS) and integrated autocorrelation time (iat) are
% computed using the initial monotone positive sequence estimator (IMSE).

if(acf_plots)
    % Compute ess and iat for all ABC-MCMC algorithms
    
    % maxlag is either 1000 or the smallest cold chain sample size
    if(single_on_multi)
        c_m = min([min(b_m(1, :)) min(b_em(1, :)) min(b_am(1, :))])-1;
    else
        c_m = min([b_em(1) b_am(1)])-1;
    end
    
    maxlag = 1000;%min([1000 c_m]);
    clear('c_m');
    imse=1; %apply IMSE when computing ESS and IAT
    maxc=6;
    
    % Standard ABC
    if run_vanilla
        [ess.vanilla, iat.vanilla] = ess_iat_maxlags(S, maxlag, imse);%, 1);
    end    
    
    if(multi_twice)
        % ABC-PTMC-K
        [ess.stdexchangem, iat.stdexchangem] = ess_iat_maxlags(S_em, maxlag, imse);%, 1);
        % ABC-APTMC-K
        [ess.anytimem, iat.anytimem] = ess_iat_maxlags(A_em, maxlag, imse);%, 1);
        if run_one
            % ABC-PTMC-1
            [ess.stdexchange, iat.stdexchange] = ess_iat_maxlags(S_e, maxlag, imse);%, 1);
            % ABC-APTMC-1
            [ess.anytime, iat.anytime] = ess_iat_maxlags(A_e, maxlag, imse);%, 1);
        end
    else
        % ABC-APTMC-K
        [ess.stdexchangem, iat.stdexchangem] = autocorr_new_multi(S_em, maxc);
        
        % ABC-APTMC-K
        [ess.anytimem, iat.anytimem] = autocorr_new_multi(A_em, maxc);
        if run_one
            % ABC-PTMC-1
            [ess.stdexchange, iat.stdexchange] = autocorr_new_multi(S_e, maxc);
            % ABC-APTMC-1
            [ess.anytime, iat.anytime] = autocorr_new_multi(A_e, maxc);
        end
    end
    
   

% COMPARE IAT %
fprintf('ABC-APTMC-1 has IAT: \n') ;
fprintf(' %.3f |', mean(mean(iat.anytime)./mean(iat.stdexchange)));
fprintf(' larger on average than ABC-PTMC-1 \n');

fprintf('ABC-PTMC-K has IAT: \n') ;
fprintf(' %.3f |', mean(mean(iat.stdexchange)./mean(iat.stdexchangem)));
fprintf(' smaller on average than ABC-PTMC-1 \n');
fprintf(' %.3f |', mean(mean(iat.anytime)./mean(iat.stdexchangem)));
fprintf(' smaller on average than ABC-APTMC-1 \n');

fprintf('ABC-APTMC-K has IAT: \n') ;
fprintf(' %.3f |', mean(mean(iat.stdexchange)./mean(iat.anytimem)));
fprintf(' smaller on average than ABC-PTMC-1 \n');
fprintf(' %.3f |', mean(mean(iat.anytime)./mean(iat.anytimem)));
fprintf(' smaller on average than ABC-APTMC-1 \n');
fprintf(' %.3f |', mean(mean(iat.stdexchangem)./mean(iat.anytimem)));
fprintf(' smaller on average than ABC-PTMC-K \n');

% COMPARE ESS % 
fprintf('ABC-APTMC-1 has ESS: \n') ;
fprintf(' %.3f |', mean(sum(ess.anytime)./sum(ess.stdexchange)));
fprintf(' larger on average than ABC-PTMC-1 \n');

fprintf('ABC-PTMC-K has ESS: \n') ;
fprintf(' %.3f |', mean(sum(ess.stdexchangem)./sum(ess.stdexchange)));
fprintf(' larger on average than ABC-PTMC-1 \n');
fprintf(' %.3f |', mean(sum(ess.stdexchangem)./sum(ess.anytime)));
fprintf(' larger on average than ABC-APTMC-1 \n');

fprintf('ABC-APTMC-K has ESS: \n') ;
fprintf(' %.3f |', mean(sum(ess.anytimem)./sum(ess.stdexchange)));
fprintf(' larger on average than ABC-PTMC-1 \n');
fprintf(' %.3f |', mean(sum(ess.anytimem)./sum(ess.anytime)));
fprintf(' larger on average than ABC-APTMC-1 \n');
fprintf(' %.3f |', mean(sum(ess.anytimem)./sum(ess.stdexchangem)));
fprintf(' larger on average than ABC-PTMC-K \n');

[sum(ess.stdexchange,1); sum(ess.anytime,1); sum(ess.stdexchangem,1); sum(ess.anytimem,1)]'
[mean(iat.stdexchange,1); mean(iat.anytime,1); mean(iat.stdexchangem,1); mean(iat.anytimem,1)]'

[mean(iat.stdexchange,1); sum(ess.stdexchange,1);-2*ones(1,3); mean(iat.anytime,1);sum(ess.anytime,1);-2*ones(1,3); 
    mean(iat.stdexchangem,1);sum(ess.stdexchangem,1);-2*ones(1,3);   mean(iat.anytimem,1); sum(ess.anytimem,1);-1*ones(1,3)]'


if all_chains_ess
    [essa.stdexchangem, iata.stdexchangem] = ess_iat_allchains(Theta_em, ne_em, n_em, imse, maxlag);
    [essa.anytimem, iata.anytimem]  = ess_iat_allchains(Theta_aem, ane_em, an_em, imse, maxlag);
    if run_one
        [essa.stdexchange, iata.stdexchange] = ess_iat_allchains(Theta_e, ne_e, n_e, imse, maxlag);
        [essa.anytime, iata.anytime] = ess_iat_allchains(Theta_ae, ane_e, an_e, imse, maxlag);
    end
    
    ESSA = zeros(params.K, 4, 3); IATA = zeros(params.K, 4, 3);
    for i=1:3
        ESSA(:,:,i) = [sum(essa.stdexchange(:,:,i),2), sum(essa.anytime(:,:,i),2),...
            sum(essa.stdexchangem(:,:,i),2), sum(essa.anytimem(:,:,i),2)];
        IATA(:,:,i) = [mean(iata.stdexchange(:,:,i),2), mean(iata.anytime(:,:,i),2),...
            mean(iata.stdexchangem(:,:,i),2), mean(iata.anytimem(:,:,i),2)];
    end
end

acf_params.maxlag = 300;
acf_params.names = {'ABC-PTMC-1', 'ABC-APTMC-1', 'ABC-PTMC-W', 'ABC-APTMC-W'};%{'Standard ABC', 'ABC-PTMC', 'ABC-APTMC'};
acf_params.linestyle = {':','--',':', '--'};
acf_params.marker = {'diamond', 'square', 'diamond', 'square'};
acf_params.markersize = {4, 3, 2, 2};
acf_params.col = {rgb('DarkBlue'), rgb('DarkOrange'),  rgb('Teal'), rgb('Red')};
[R_out]=acf_keep_maxlag({S_e, A_e, S_em, A_em}, {b_e(1,:),b_a(1, :), b_em(1,:),b_am(1, :)}, acf_params.maxlag);
figure;
subplot(3, 1, 1); plot_acf_multi(R_out, 1, acf_params, '\theta_1', 1);
subplot(3, 1, 2); plot_acf_multi(R_out, 2, acf_params, '\theta_2', 1);
subplot(3, 1, 3); plot_acf_multi(R_out, 3, acf_params, '\theta_3', 1);
end

%% Chain plots
% Plotting the chains output by the various ABC-MCMC algorithms

if(chain_plots)
    % Multi processor chains
    if(run_vanilla)
        plot_multi_chains(S, true_theta, 'Standard ABC on K processors')
    end
    if(multi_twice)
        plot_multi_chains(S_em, true_theta, 'ABC-PTMC-K')
        plot_multi_chains(A_em, true_theta, 'ABC-APTMC-K')
        if run_one
            plot_multi_chains(S_e, true_theta, 'ABC-PTMC-1')
            plot_multi_chains(A_e, true_theta, 'ABC-APTMC-1')
        end
    else
        plot_chains(S_em, true_theta, 'ABC-PTMC-K')
        plot_chains(A_em, true_theta, 'ABC-APTMC-K')
    end
end

%% Marginal plots
% Plot histograms of the marginal posteriors output from ABC-MCMC
% algorithms and compare them to the ABC rejection sampling reference

if(marginal_plots)
    nbins=50;
    if run_vanilla
        plot_marginal_posteriors(cell2mat(S'), S_rej, nbins, 'Standard ABC on K processors');
    end
    if(multi_twice)
        plot_marginal_posteriors(cell2mat(S_em'), S_rej, nbins, 'ABC-PTMC-K');
        plot_marginal_posteriors(cell2mat(A_em'), S_rej, nbins, 'ABC-APTMC-K');
        if run_one
            plot_marginal_posteriors(cell2mat(S_e'), S_rej, nbins, 'ABC-PTMC-1');
            plot_marginal_posteriors(cell2mat(A_e'), S_rej, nbins, 'ABC-APTMC-1');
        end
    else
        %only adapted for 4 processors
        plot_marginal_posteriors(S_em, S_rej, nbins, 'ABC-PTMC-K');
        plot_marginal_posteriors(A_em, S_rej, nbins, 'ABC-APTMC-K');
    end
end

if(density_plots&&run_anytime)
    % Display all marginal densities on the same plot
    den = [2^8 2^16 2^8
        2^7 2^7 2^7
        2^7 2^7 2^7
        2^7 2^7 2^7];
    
    if(multi_twice)
        plot_three_densities(cell2mat(S_e'),...
            cell2mat(S_em'),...
            cell2mat(A_em'),...
            S_rej, den);
    else
        plot_three_densities([S{1}; S{2}; S{3}; S{4}], S_em, A_em, S_rej, den);
    end
end

if(mean_plots&&run_anytime)
    % Display the evolution of the sample means for all parameters
    plot_three_means(A_e{1}, S_em{1}, A_em{1}, S_rej);
end

%% Timeplots (ABC-PTMC-1)
% Timelines of local and exchange moves for the ABC-PTMC algorithm for the
% first 500 seconds of their runs, to observe the time the algorithm spends
% performing each type of update.

if(time_plots)
    % Timeline plots (one processor)
    time_interval = [0 min(300, params.T)]; %[params.T-500 params.T];
    
    % Timeline plots (K processors)
    names= {'Global', 'Worker 1', 'Worker 2', 'Worker 3','Worker 4'};
    %{'Updates within workers', 'Exchange between workers',...
    %    'Worker 1: local', 'Worker 1: exchange', 'Worker 2: local', 'Worker 2: exchange',...
    %    'Worker 3: local', 'Worker 3: exchange', 'Worker 4: local', 'Worker 4: exchange'};
    col = {rgb('MediumBlue'), rgb('SkyBlue'), rgb('SkyBlue'),...
        rgb('SkyBlue'), rgb('SkyBlue')};
    edgecol = {rgb('SkyBlue'), rgb('SteelBlue'), rgb('SteelBlue'), rgb('SteelBlue'), rgb('SteelBlue'),...
        rgb('SteelBlue'), rgb('SteelBlue'), rgb('SteelBlue'), rgb('SteelBlue')};
    multi_timeplots_twice(TM_m{1}, 'ABC-PTMC-W: Timeline of local and exchange moves', time_interval, names,1 , col, edgecol);
    multi_timeplots_twice(TM_am{1}, 'ABC-APTMC-W: Timeline of local and exchange moves', time_interval, names, 0);
    
    
    % Histograms (for an overview of all local/exchange move times)
    % figure; histogram(TM(TM(:,2)==1,1), 'NumBins', 300,'DisplayName', 'ABC-PTMC-1 local move times');
    % hold on; xlabel('local move times');
    % histogram(TM(TM(:,2)==2,1), 'NumBins', 80,'DisplayName', 'ABC-PTMC-1 exchange move times');
end

%% Boxplots (ABC-PTMC-1 and ABC-APTMC-1)
% Boxplots of all the times taken to complete local moves (since displaying
% them all on a timeline would make it unreadable). Requires anytime to
% have run

if(box_plots&&run_anytime)
    % Boxplots of local move times
    
    if(single_on_multi)
        w=1; TMS = TM{w};
        if(run_anytime)
            TMA = TM_a{w};
        end
    else
        TMS = TM;
        if(run_anytime)
            TMA = TM_a;
        end
    end
    
    % Create ‘NaN’ Array
    C3 = NaN(max([sum(TMS(:,2)==1), sum(TMA(:,2)==1)]), 2);
    C3(1:sum(TMS(:,2)==1),1) = TMS(TMS(:,2)==1,1);
    C3(1:sum(TMA(:,2)==1), 2) =  TMA(TMA(:,2)==1,1);
    figure; bxplt=boxplot(C3, {'ABC-PTMC-1', 'ABC-APTMC-1'});
    ylabel('local move times (seconds)');
    bxplt.YLim = [0 10];
    title('Boxplot of time taken to perform local moves');
end
% 
fprintf('Saving workspace... ')
  save(sprintf('results/ABC/LV/LV_example_wkspace_multi_%s.mat', datestr(now, 'yyyymmdd_HHMM')))
fprintf('done \n')
% 
% save (or read) everything
% make folder
% algorithm_names = {'anytime', 'stdexchange', 'stdexchangem','anytimem'};
% saving = 0;
% for nalgs=1:length(algorithm_names)
%     algorithm_name = algorithm_names{nalgs};
%     dirname =  sprintf('LV_multi_20181203_%s_%d', algorithm_name, params.T);
%     %     create folder to store results
%     if(saving)
%         mkdir(sprintf('results/ABC/LV/%s', dirname))
%     end
%     switch algorithm_name
%         case 'anytime'
%             if(saving)
%                 save_all(dirname, algorithm_name, A_e, b_a, TM_a, out_a, rates_a)
%             else
%                 %[A_e, b_a, TM_a, out_a, rates_a] = read_all(dirname, algorithm_name, A_e, b_a, TM_a, out_a, rates_a);
%                 [A_e, b_a, TM_a, out_a, rates_a] = read_all(dirname, algorithm_name, ntimes);                
%             end
%             
%         case 'anytimem'
%             if(saving)
%                 save_all(dirname, algorithm_name, A_em, b_am, TM_am, out_am, rates_am)
%             else
%                 %[A_em, b_am, TM_am, out_am, rates_am] = read_all(dirname, algorithm_name, A_em, b_am, TM_am, out_am, rates_am);
%                 [A_em, b_am, TM_am, out_am, rates_am] = read_all(dirname, algorithm_name, ntimes);
% 
%             end
%             
%         case 'stdexchange'
%             if(saving)
%                 save_all(dirname, algorithm_name, S_e, b_e, TM, out_e, rates_e)
%             else
%                 %[S_e, b_e, TM, out_e, rates_e] = read_all(dirname, algorithm_name, S_e, b_e, TM, out_e, rates_e);
%                 [S_e, b_e, TM, out_e, rates_e] = read_all(dirname, algorithm_name, ntimes);
%             end
%             
%         case 'stdexchangem'
%             if(saving)
%                 save_all(dirname, algorithm_name, S_em, b_em, TM_m, out_em, rates_em)
%             else
%                 %[S_em, b_em, TM_m, out_em, rates_em] = read_all(dirname, algorithm_name, S_em, b_em, TM_m, out_em, rates_em);
%                 [S_em, b_em, TM_m, out_em, rates_em] = read_all(dirname, algorithm_name, ntimes);
%             end
%             
%         otherwise
%             fprintf('No algorithm with that name available')
%     end
%     
% end

% delete(gcp);
oneto20=1:20;
tbl = [oneto20', zeros(20,1)-2, round(params.epsilon,3)',zeros(20,1)-2, round(sigma, 3)',zeros(20,1)-2, mean(b_e,2),zeros(20,1)-2, mean(b_em,2), zeros(20,1)-2, mean(b_a,2), zeros(20,1)-2, mean(b_am,2), -1*ones(20,1)]
tbl2 = [oneto20', zeros(20,1)-2, round(params.epsilon,3)',zeros(20,1)-2, round(sigma, 3)',zeros(20,1)-2, mean(sum(out_e(:,4:5,:),2),3),zeros(20,1)-2, mean(sum(out_em(:,4:5,:),2),3), zeros(20,1)-2, mean(sum(out_a(:,4:5,:),2),3), zeros(20,1)-2, mean(sum(out_am(:,4:5,:),2),3), -1*ones(20,1)]

tbl = [oneto20', round(params.epsilon,3)', round(sigma, 3)', mean(b_e,2), mean(b_em,2), mean(b_a,2), mean(b_am,2)]
tbl2 = [oneto20', round(params.epsilon,3)', round(sigma, 3)', mean(sum(out_e(:,4:5,:),2),3), mean(sum(out_em(:,4:5,:),2),3), mean(sum(out_a(:,4:5,:),2),3), mean(sum(out_am(:,4:5,:),2),3)]

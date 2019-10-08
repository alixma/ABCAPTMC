%% Initialisation
% Applying approximate Bayesian computation (ABC) with parallel tempering
% Monte Carlo (PTMC) and anytime parallel tempering Monte Carlo (APTMC) to
% estimate the parameters $\theta = \left(\theta_1, \theta_2,
% \theta_3\right)$ of a stochastic Lotka-Volterra predator-prey model.

% Various flags to decide which algorithms and diagnostics are run
%clear;
run_rejection = 0; run_vanilla =1; run_stdexchange =1; run_anytime=1;
run_one = 1; single_on_multi = 1; run_multi=0; oldess=0;
multi_twice = 0; essovertime=1;
% Plots and diagnostics
acf_plots = 1; time_plots = 0; box_plots = 0;
density_plots = 0; mean_plots = 0;
chain_plots = 0; marginal_plots = 0;


LV_initialisation
fprintf('I will be done on:')
disp(datetime(clock + [0 0 0 0 0 (run_vanilla+run_anytime+run_stdexchange)*params.T]))

%% Standard ABC-MCMC algorithm
% The vanilla MCMC algorithm which employs ABC in its local moves. The
% algorithm runs for on a chain associated with the ball of
% radius $\varepsilon = 1$.

if(run_vanilla)
    % Run LV_standard only performing Single processor ABC-MCMC
    standard_rejection=0; 
    run_single_standard=1-single_on_multi; 
    run_multi_standard=single_on_multi;
    LV_standard    
end

if(run_one)
    %% Standard ABC-MCMC algorithm with exchange moves (ABC-PTMC-1 and ABC-PTMC-K)
    % The above algorithm but run sequentially on $K$ chains, each targetting a posterior
    % associated with a different $\varepsilon^k$. The first, _cold_ chain is the
    % one associated with ball of radius $\varepsilon^1 = 1$, the _warmer_
    % chains target posteriors associated with balls of larger radii. Every
    % $\delta_t$ local moves an exchange move occurs between two uniformly
    % randomly chosen adjacent chains.
    
    if(run_stdexchange)
        % Run LV_stdexchange only performing single processor ABC-PTMC-1
        run_single_stdexchange=1-single_on_multi;
        run_single_stdexchange_on_multi=single_on_multi;
        run_multi_stdexchange=0; run_multi_stdexchange_twice=0;
        LV_stdexchange
    end
    
    %% Anytime ABC-MCMC algorithm with exchange moves (ABC-APTMC-1 and ABC-APTMC-K)
    % The anytime version of the ABC_APTMC-1 algorithm, corrected for length bias.
    % Again run on $K=8$ chains, each targetting a posterior associated with a
    % different $\varepsilon^k$. This time, exchange moves occur every
    % $\delta^a_t$ seconds, where now $\delta^a_t$ is the median time the
    % previous algorithm took to perform $\delta_t$ local moves.
    
    if(run_anytime)
        % Run LV_anytime only performing single processor ABC-APTMC-1
        run_single_anytime=1-single_on_multi;
        run_single_anytime_on_multi=single_on_multi;
        run_multi_anytime=0; run_multi_anytime_twice=0;
        LV_anytime
    end
end

if(run_multi)
    
   
    
    %% Multi-processor standard ABC-MCMC algorithm with exchange moves (ABC-PTMC-W)
    % The above algorithm but run in parallel on 4 workers, with $Kk=2$
    % chains per processor (i.e. K=8 chains in total), each chain targets a
    % posterior associated with a different $\varepsilon^k$ and exchange moves
    % are performed between workers.
    
    if(run_stdexchange)
        
        % Run LV_stdexchange only performing multi-processor ABC-PTMC-K        
        run_single_stdexchange=0; run_single_stdexchange_on_multi=0;
        run_multi_stdexchange=1-multi_twice; run_multi_stdexchange_twice=multi_twice;
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
        run_multi_anytime=1-multi_twice; run_multi_anytime_twice=multi_twice;
        LV_anytime
        
    end
end


%% Performance comparison
% All ABC-MCMC algorithms ran for the same amount of time. To compare their
% efficiency, sample autocorelation function plots are drawn and the
% effective sample size (ESS) and integrated autocorrelation time (iat) are
% computed using the initial monotone positive sequence estimator (IMSE).

if(acf_plots)
    % Compute ess and iat for all ABC-MCMC algorithms
    
    % maxlag is either 1000 or the smallest cold chain sample size
    if(run_one)
        if(single_on_multi)
            c = min([min(b_m(1, :)) min(b_e(1, :)) min(b_a(1, :))])-1;
        else
            c = min([b b_e(1) b_a(1)])-1;
        end
    else
        c = 1000;
    end
    if(run_multi)
        if(single_on_multi)
            c_m = min([min(b_m(1, :)) min(b_em(1, :)) min(b_am(1, :))])-1;
        else
            c_m = min([b_em(1) b_am(1)])-1;
        end
    else
        c_m = 1000;
    end

    maxlag = min([1000 c c_m]);
    clear('c', 'c_m');
    imse=1; %apply IMSE when computing ESS and IAT
    maxc = 6;
    
    if(run_one)
        if(oldess)
            vanilla_ess = zeros(W, 3); vanilla_iat = zeros(W, 3);
            stdexchange_ess = zeros(W, 3); stdexchange_iat = zeros(W, 3);
            anytime_ess = zeros(W, 3); anytime_iat = zeros(W, 3);
            for w=2:W
                % Standard ABC
                S_w = S{w};
                [vanilla_ess(w, 1), vanilla_iat(w, 1)] = ESS_IAT(S_w(:, 1), maxlag, imse);
                [vanilla_ess(w, 2), vanilla_iat(w, 2)] = ESS_IAT(S_w(:, 2), maxlag, imse);
                [vanilla_ess(w, 3), vanilla_iat(w, 3)] = ESS_IAT(S_w(:, 3), maxlag, imse);
                
                % ABC-PTMC-1
                S_e_w = S_e{w};
                [stdexchange_ess(w, 1), stdexchange_iat(w, 1)] = ESS_IAT(S_e_w(:, 1), maxlag, imse);
                [stdexchange_ess(w, 2), stdexchange_iat(w, 2)] = ESS_IAT(S_e_w(:, 2), maxlag, imse);
                [stdexchange_ess(w, 3), stdexchange_iat(w, 3)] = ESS_IAT(S_e_w(:, 3), maxlag, imse);
                
                % ABC-APTMC-1
                A_e_w = A_e{w};
                [anytime_ess(w, 1), anytime_iat(w, 1)] = ESS_IAT(A_e_w(:, 1), maxlag, imse);
                [anytime_ess(w, 2), anytime_iat(w, 2)] = ESS_IAT(A_e_w(:, 2), maxlag, imse);
                [anytime_ess(w, 3), anytime_iat(w, 3)] = ESS_IAT(A_e_w(:, 3), maxlag, imse);
            end
            
            % Standard ABC
            ess.standard1 = sum(vanilla_ess(:,1)); iat.standard1 = mean(vanilla_iat(:,1));
            ess.standard2 = sum(vanilla_ess(:,2)); iat.standard2 = mean(vanilla_iat(:,2));
            ess.standard3 = sum(vanilla_ess(:,3)); iat.standard3 = mean(vanilla_iat(:,3));
            
            % ABC-PTMC-1
            ess.stdexchange1 = sum(stdexchange_ess(:,1)); iat.stdexchange1 = mean(stdexchange_iat(:,1));
            ess.stdexchange2 = sum(stdexchange_ess(:,2)); iat.stdexchange2 = mean(stdexchange_iat(:,2));
            ess.stdexchange3 = sum(stdexchange_ess(:,3)); iat.stdexchange3 = mean(stdexchange_iat(:,3));
            
            % ABC-APTMC-1
            ess.anytime1 = sum(anytime_ess(:,1)); iat.anytime1 = mean(anytime_iat(:,1));
            ess.anytime2 = sum(anytime_ess(:,2)); iat.anytime2 = mean(anytime_iat(:,2));
            ess.anytime3 = sum(anytime_ess(:,3)); iat.anytime3 = mean(anytime_iat(:,3));
            clear('S_m_w', 'S_e_w', 'A_e_w');
        else
            % Standard ABC
            [ess.standard1, iat.standard1] = autocorr_new(S, 1, maxc);
            [ess.standard2, iat.standard2] = autocorr_new(S, 2, maxc);
            [ess.standard3, iat.standard3] = autocorr_new(S, 3, maxc);
            % ABC-PTMC-1
            [ess.stdexchange1, iat.stdexchange1] = autocorr_new(S_e, 1, maxc);
            [ess.stdexchange2, iat.stdexchange2] = autocorr_new(S_e, 2, maxc);
            [ess.stdexchange3, iat.stdexchange3] = autocorr_new(S_e, 3, maxc);
            
            if(run_anytime)
                % ABC-APTMC-1
                [ess.anytime1, iat.anytime1] = autocorr_new(A_e, 1, maxc);
                [ess.anytime2, iat.anytime2] = autocorr_new(A_e, 2, maxc);
                [ess.anytime3, iat.anytime3] = autocorr_new(A_e, 3, maxc);
                
            end
        end
    end    
    
    ess
    iat
    
    standard.ess = [ess.standard1 ess.standard2 ess.standard3];
    standard.iat = [iat.standard1 iat.standard2 iat.standard3];
    
    if(run_one)
        stdexchange.ess = [ess.stdexchange1 ess.stdexchange2 ess.stdexchange3];
        stdexchange.iat = [iat.stdexchange1 iat.stdexchange2 iat.stdexchange3];
        anytime.ess = [ess.anytime1 ess.anytime2 ess.anytime3];
        anytime.iat = [iat.anytime1 iat.anytime2 iat.anytime3];
        
        
        fprintf('Standard ABC has IAT: \n') ;
        fprintf(' %.3f larger on average than ABC-PTMC-1 \n', mean(standard.iat./stdexchange.iat));
        fprintf(' %.3f larger on average than ABC-APTMC-1 \n', mean(standard.iat./anytime.iat));
        
        fprintf('Standard ABC has ESS: \n') ;
        fprintf(' %.3f smaller on average than ABC-PTMC-1 \n', mean(stdexchange.ess./standard.ess));
        fprintf(' %.3f smaller on average than ABC-APTMC-1 \n', mean(anytime.ess./standard.ess));
    end
   
    
    if(run_one)
        % Plot the sample acf for the ABC, ABC-PTMC-1 and ABT-APTMC-1 algorithms        
%         figure; hold on;
%         names = {'Standard ABC', 'ABC-PTMC-1', 'ABC-APTMC-1'};
%         subplot(3, 1, 1);plot_three_acf(S, S_e, A_e, 1, '\theta_1', names, maxlag, single_on_multi);
%         subplot(3, 1, 2);plot_three_acf(S, S_e, A_e, 2, '\theta_2', names, maxlag, single_on_multi);
%         subplot(3, 1, 3);plot_three_acf(S, S_e, A_e, 3, '\theta_3', names, maxlag, single_on_multi);
        
        acf_params.maxlag = 300;
        acf_params.names = {'Standard ABC', 'ABC-PTMC', 'ABC-APTMC'};
        acf_params.linestyle = {':','--',':', '--'};
        acf_params.marker = {'diamond', 'square', 'diamond', 'square'};
        acf_params.markersize = {4, 3, 2, 2};
        acf_params.col = { [12 195 82] ./ 255, rgb('DarkBlue'),  rgb('DarkOrange'), rgb('Red')};
        [R_out]=acf_keep_maxlag({S, S_e, A_e}, {b_m(1,:), b_e(1,:),b_a(1, :)}, acf_params.maxlag);
        figure;
        subplot(3, 1, 1); plot_acf_multi(R_out, 1, acf_params, '\theta_1', 1);
        subplot(3, 1, 2); plot_acf_multi(R_out, 2, acf_params, '\theta_2', 1);
        subplot(3, 1, 3); plot_acf_multi(R_out, 3, acf_params, '\theta_3', 1);
    end
   
end

%% Chain plots
% Plotting the chains output by the various ABC-MCMC algorithms

if(chain_plots)
    % Single processor chains
    if(run_one)
        if(single_on_multi)
            plot_multi_chains(S, true_theta, 'Standard ABC on K processors')
            plot_multi_chains(S_e, true_theta, 'ABC-PTMC-1 on K processors')
            plot_multi_chains(A_e, true_theta, 'ABC-APTMC-1 on K processors')
        else
            plot_chains(S, true_theta, 'Standard ABC')
            plot_chains(S_e, true_theta, 'ABC-PTMC-1')
            plot_chains(A_e, true_theta, 'ABC-APTMC-1')
        end
    end
    
end

%% Marginal plots
% Plot histograms of the marginal posteriors output from ABC-MCMC
% algorithms and compare them to the ABC rejection sampling reference

if(marginal_plots)
    nbins=70;
    if(run_one)
        if(single_on_multi)
            %only adapted for 4 processors
            plot_marginal_posteriors(cell2mat(S'), S_rej, nbins, 'Standard ABC on K processors')
            plot_marginal_posteriors(cell2mat(S_e'), S_rej, nbins, 'ABC-PTMC-1 on K processors')
            plot_marginal_posteriors(cell2mat(A_e'), S_rej, nbins, 'ABC-APTMC-1 on K processors')
        else
            plot_marginal_posteriors(S, S_rej, nbins, 'Standard ABC');
            plot_marginal_posteriors(S_e, S_rej, nbins, 'ABC-PTMC-1');
            plot_marginal_posteriors(A_e, S_rej, nbins, 'ABC-APTMC-1');
        end
    end    
end

if(density_plots&&run_anytime)
    % Display all marginal densities on the same plot
    den = [2^8 2^16 2^8
        2^8 2^8 2^8
        2^8 2^8 2^8
        2^8 2^8 2^8];
    if(single_on_multi)
        plot_three_densities(cell2mat(S_e'),...
            cell2mat(A_e'),...
            cell2mat(S'),...
            S_rej, den);
    else
        plot_three_densities(S, S_e, A_e, S_rej, den);
    end
    
end

if(mean_plots&&run_anytime)
    w=2;
    % Display the evolution of the sample means for all parameters
    plot_three_means(S{2}, S_e{2}, A_e{2}, S_rej);
end

%% Timeplots (ABC-PTMC-1)
% Timelines of local and exchange moves for the ABC-PTMC algorithm for the
% first 500 seconds of their runs, to observe the time the algorithm spends
% performing each type of update.

if(time_plots)
    % Timeline plots (one processor)
    time_interval = [0 min(300, params.T)]; %[params.T-500 params.T];
    if(run_one)
        if(single_on_multi)
            w=4; TMS = TM{w};
            if(run_anytime)
                TMA = TM_a{w};
            end
        else
            TMS = TM;
            if(run_anytime)
                TMA = TM_a;
            end
        end
        figure;
        subplot(2, 1, 1);
        timeplots(TMS(1:300, :), 'ABC-PTMC-1: Timeline of local and exchange moves', {rgb('MediumBlue'), rgb('LightCyan')}, {rgb('LightCyan'), rgb('MediumBlue')}, time_interval);
        subplot(2, 1, 2);
        timeplots(TMA(1:300, :), 'ABC-APTMC-1: Timeline of local and exchange moves', {rgb('DarkOrange'), rgb('SandyBrown')}, {rgb('DarkRed'), rgb('DarkOrange')}, time_interval);
    end  
    
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

if essovertime
     %w=5; 
    TMA = []; TMS = []; TMM=[];
    for w=1:W
        TMA = [TMA; TM_a{w}(:, 3:4)]; 
        TMS = [TMS; TM{w}(:, 3:4)]; 
        TMM = [TMM; TM_m{w}(TM_m{w}(:, 1)>0, :)]; 
    end
    
    iat_a = [iat.anytime1, iat.anytime2, iat.anytime3];%[28.9083 52.6397 71.7799];
    iat_e = [iat.stdexchange1, iat.stdexchange2, iat.stdexchange3];%[12.4816 14.5992 17.2530];
    iat_m = [iat.standard1, iat.standard2, iat.standard3];%[49.4429 98.1508 102.9551];
    col = [[12 195 82] ./ 255; rgb('DarkBlue');  rgb('Red')];
    scol = [rgb('DarkSeaGreen'); rgb('LightBlue');  rgb('Goldenrod')];
    figure;
    sz=.1; alph=.5; lw=1.5;
    for i=1:3
        p_a = fitlm(TMA(:, 1), TMA(:, 2)/iat_a(i), 'Intercept', false); [fit_a, int_a] = predict(p_a, TMA(:, 1));
        p_e = fitlm(TMS(:, 1), TMS(:, 2)/iat_e(i), 'Intercept', false); [fit_e, int_e] = predict(p_e, TMS(:, 1));
        p_m = fitlm(TMM(:, 1), TMM(:, 2)/iat_m(i), 'Intercept', false); [fit_m, int_m] = predict(p_m, TMM(:, 1));
        subplot(3, 1, i);
        hold on;
        g1 = scatter(TMA(:, 1), TMA(:, 2)/iat_a(i), sz, 'Marker', '.',...
                                    'MarkerEdgeColor', scol(3,:),...
                                    'MarkerFaceColor', scol(3,:),...
                                    'MarkerFaceAlpha', alph/5,...
                                    'MarkerEdgeAlpha', alph/5,...
                                    'DisplayName','ABC-APTMC');
        g2 = scatter(TMS(:, 1), TMS(:, 2)/iat_e(i), sz, 'MarkerEdgeColor', scol(2,:),...
                                    'MarkerFaceColor', scol(2,:),...
                                    'MarkerFaceAlpha', alph,...
                                    'MarkerEdgeAlpha', alph,... 
                                    'DisplayName','ABC-PTMC');
        g3 = scatter(TMM(:, 1), TMM(:, 2)/iat_m(i), sz, 'Marker', 'd',...
                                    'MarkerEdgeColor', scol(1,:),...
                                    'MarkerFaceColor', scol(1,:),...
                                    'MarkerFaceAlpha', alph,...
                                    'MarkerEdgeAlpha', alph,...
                                    'DisplayName','Standard');             
        h1 = plot(TMA(:, 1), fit_a, 'Color', col(3,:), 'LineWidth', lw, 'DisplayName','ABC-APTMC fit');
        %h1_1 = plot(TMA(:, 1), int_a(:, 1), 'LineStyle', '--', 'Color', col(3,:), 'DisplayName','ABC-APTMC CI');
        %h1_2 = plot(TMA(:, 1), int_a(:, 2), 'LineStyle', '--', 'Color', col(3,:));
        h2 = plot(TMS(:, 1), fit_e, 'Color', col(2,:), 'LineWidth', lw, 'DisplayName','ABC-PTMC fit');
        %h2_1 = plot(TMS(:, 1), int_e(:, 1), 'LineStyle', '--', 'Color', col(2,:), 'DisplayName','ABC-PTMC CI');
        %h2_2 = plot(TMS(:, 1), int_e(:, 2), 'LineStyle', '--', 'Color', col(2,:));
        h3 = plot(TMM(:, 1), fit_m, 'Color', col(1,:), 'LineWidth', lw, 'DisplayName','Standard fit');
        %h3_1 = plot(TMM(:, 1), int_m(:, 1), 'LineStyle', '--', 'Color', col(1,:), 'DisplayName','Standard CI');
        %h3_2 = plot(TMM(:, 1), int_m(:, 2), 'LineStyle', '--', 'Color', col(1,:));
        xlabel('time'); xlim([0 T_s]);
        ylabel('ESS')
        title(sprintf('ESS over time, \\theta_%d', i))
        legend('Location', 'northwest');%([h1 h1_1 h2 h2_1 h3 h3_1])
    end
end

fprintf('Saving workspace... ')
save(sprintf('results/ABC/LV/LV_example_wkspace_%s.mat', datestr(now, 'yyyymmdd_HHMM')))
fprintf('done \n')

delete(gcp);
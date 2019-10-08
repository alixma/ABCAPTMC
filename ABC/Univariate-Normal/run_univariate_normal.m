%% ABC:UNIVARIATE NORMAL EXAMPLE %%
observations.y = 3;
observations.sigma = 1;
params.prior = [0 sqrt(5)];
params.posterior = [5/2 sqrt(5/6)];
params.rho = 0.5; pre=1; 
exchange.simple=1; 
if exchange.simple
    exmoves = "simple";
else
    exmoves = "";
end
params.K=10; Lambda = params.K;
timepar.burnin=30; timepar.T = timepar.burnin+3600;  
params.N=5e6; exchange.deltat=5e-4;
params.epsilon = logspace(log10(0.1), log10(1.1), params.K);

C = 4;
if(isempty(gcp('nocreate')))
     parpool(min(C, 4));
end

%% Run ABC rejection samplers
fprintf("Simulating a few samples from posterior...")
[S_rej, n_rej] = ABC_rejection(1, params, observations);
fprintf("done\n")

%% Run ABC-APTMC algortithm C times
fprintf("ABC-APTMC algorithm run %d times on %d chains for %d seconds\n", C, params.K, timepar.T)
fprintf("  %s exchange moves every %.4f seconds\n", exmoves, exchange.deltat)
params.theta_in = reshape(randsample(S_rej, C*params.K, 0), params.K, C);
% corrected
fprintf("  corrected for bias\n")
exchange.correct=1;
tic
[Theta, X, Rej, n, ne, nsw, sw, TM] = ABC_PT_AMC(C, timepar, params, exchange, observations);
toc

% Get rid of samples from burnin
[R, b] = remove_burnin(Rej, ne, n, 0);
c=1;
[out, rates] = summary_chain({R{c}, b, nsw, sw}, {params.epsilon, params.rho}, 0);

[S, ~] = remove_burnin(Theta, ne, n, 0);

% Plot all chains
plot_allchains(S, params, observations, 40, 2^8)
fprintf('Chains composition (3->7) and rates (8->9): \n')
disp([(1:params.K)' round(params.epsilon', 2) round(out, 3) round(rates, 3)])

fprintf("\n___________________________________________________________\n")

fprintf("ABC-APTMC algorithm run %d times on %d chains for %d seconds\n", C, params.K, timepar.T)
fprintf("  %s exchange moves every %.4f seconds\n", exmoves, exchange.deltat)
%% uncorrected
fprintf("  uncorrected for bias\n")
exchange.correct = 0;
tic
[Theta1, X1, Rej1, n1, ne1, nsw1, sw1, TM1] = ABC_PT_AMC(C, timepar, params, exchange, observations);
toc

% Get rid of samples from burnin
[R1, b] = remove_burnin(Rej1, ne1, n1, 0);
c=1;
[out1, rates1] = summary_chain({R1{c}, b, nsw, sw}, {params.epsilon, params.rho}, 0);
fprintf('Chains composition (3->7) and rates (8->9): \n')
disp([(1:params.K)' round(params.epsilon', 2) round(out1, 3) round(rates1, 3)])

[S1, ~] = remove_burnin(Theta1, ne1, n1, 0);

% Plot all chains
plot_allchains(S1, params, observations, 40, 2^8)

%save workspace
fprintf('Saving workspace... ')
save(sprintf('results/univariate_normal_wkspace_%s.mat', datestr(now, 'yyyymmdd_HHMM')))
fprintf('done \n')

%%
% combine parallel runs for plotting both algorithms together
if C>1
    S_u = cell(1, params.K);
    S_c = cell(1, params.K);
    
    for k=1:params.K
        S_u{k} = [];
        S_c{k} = [];
        for c=1:C
            S_u{k} = [S_u{k} S1{c}{k}'];
            S_c{k} = [S_c{k} S{c}{k}'];
        end

    end
else
    S_u = S1;
    S_c = S;
end

plot_allchains_both({S_c, S_u}, params, observations, 2^10)

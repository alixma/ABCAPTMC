%% Initialisation
% Applying approximate Bayesian computation (ABC) with parallel tempering
% Monte Carlo (PTMC) and anytime parallel tempering Monte Carlo (APTMC) to
% estimate the parameters $\theta = \left(\theta_1, \theta_2,
% \theta_3\right)$ of a stochastic Lotka-Volterra predator-prey model.

% observations
true_theta = [1 0.005 0.6];
observations.y = [88 165 274 268 114 46 32 36 53 92];
observations.ET = 10 + 1; observations.dt = 1;

% Lotka-Volterra model settings
LV.M=[50; 100];
LV.Pre = [1 0; 1 1; 0 1];
LV.Post = [2 0; 0 2; 0 0];
LV.h = @(y, th) [th(1)*y(1), th(2)*y(1)*y(2), th(3)*y(2)];

% Various parameters
% Run time of algorithm
params.burnin = 3600;
params.T = params.burnin+9*3600; params.N=1e5;
% prior
params.prior = [1 1 1]; %[1 100 1]; %[1 0.01 1];
params.pr = @(th) 1 ;

nsweeps = 3;
ntimes = 1;

%% ABC Rejection sampling
% For reference, obtain an approximation of the posterior by generating
% $N$ samples from the prior, simulating observations and only retaining
% those which hit the ball of radius $\varepsilon = 1$.

if(run_rejection)
    % Run LV_standard only performing ABC rejection sampling
    standard_rejection = 1; run_single_standard = 0; run_multi_standard = 0;
    LV_standard
else
    % If already run, fetch results from files
    S_rej = dlmread(sprintf('results/ABC/LV/LV_rejection_%d.csv', 1988));
    %S_rej = dlmread(sprintf('results/ABC/LV/LV_rejection_%d.csv', 2364));
    params.S_rej =  S_rej;
    fprintf('Rejection ABC: Acceptance rate = %f percent \n', 1988/1e8*100)
end
%%
% This algorithm has an extremely low acceptance rate and is therefore
% highly inefficient, which is why it is a good idea to resort to MCMC.

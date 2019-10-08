function [Theta_out, n_out, Theta, X, Rej, n] = ABC_rejection(C, params, observations)
% Single Processor ABC-APTMC algorithm
% Note:         Outputs all chains
% Inputs:       K: # chains
%               Tt: {global deadline, 'real' or 'virtual'}
%               N:  # slots allocated to matrix X to record samples
%               C:  run the algorithm C times in parallel
%               y: observation(s)
%               sigma: standard deviation for simulating from likelihood
%               epsilon: matrix of ball radii for all chains (tolerance levels)
%               prior_params: parameters of Gaussian prior
%               rho: standard deviation for random walk metropolis proposal
%               deltat: time between exchange moves
%               pre: include prior check (true (1) by default)
%               correct: apply bias correction (1) or not (0)
% Outputs:      Theta, X: cold chain(s) and corresponding simulations
%               n: sample size obtained for each chain 
%               ne: first sample at minimum epsilon for each chain 
%               nsw: total # exchange moves performed   
%               sw: total # exchange moves accepted   

N = params.N; E = params.epsilon(end, 1);

Theta = zeros(N, C);
Theta_out = cell(C,1);
X = zeros(N, C);
Rej = zeros(N, C);
n = ones(1, C);


for c=1:C
    par = params; obs=observations;
    % initialise theta
    Theta1 = normrnd(par.prior(1), par.prior(2), N, 1);
    X1 = normrnd(Theta1, obs.sigma);
    Rej1 = abs(obs.y - X1)<E;
    
    Theta(:,c) = Theta1;
    X(:,c) = X1;
    Rej(:,c) = Rej1;  
    n(c) = sum(Rej1);
    Theta_out{c} = Theta1(Rej1);
end
Theta_out = cell2mat(Theta_out);
n_out = sum(n);
fprintf("Acceptance rate: %.2f%% \n", n_out/N)
end
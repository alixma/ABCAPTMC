function [Theta_out, X_out, n] = LV_ABC_rejection(params, LV, observations, epsilon)
% Rejection ABC algorithm for Lotka-Volterra example
% --> used to obtain initial samples and posterior approximations
N = params.T;

% Sample from prior
expmeans = repmat(1./params.prior, N, 1);
Theta = exprnd(expmeans, [N, 3]);


% initialise theta and x for each chain
deadline.epsilon = epsilon;
deadline.T = 24*3600;
deadline.all = tic;
[X, ~] = get_x_multi(N, LV, observations, Theta, deadline, 1);

Y = zeros(N, 3);
Y(:,1) = sum(abs(log(X) - log(observations.y1))<1, 2);

keep = (Y(:,3)==10);

Theta_out = Theta(keep,:);
X_out = X(keep,:);
n = sum(keep);
end
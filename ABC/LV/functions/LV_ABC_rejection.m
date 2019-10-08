function [Theta_out, X_out, n] = LV_ABC_rejection(params, LV, observations, epsilon)
% Rejection ABC algorithm for Lotka-Volterra example
N=params.T;

% Sample from prior
if isfield(observations, 'y2')
    expmeans = repmat(params.prior, N, 1); %repmat(1./params.prior, N, 1);
    Theta = unifrnd(0, expmeans);
else
    expmeans = repmat(1./params.prior, N, 1);
    Theta = exprnd(expmeans, [N, 3]);
end

% initialise theta and x for each chain
deadline.epsilon = epsilon;
deadline.T = 24*3600;
deadline.all = tic;
[X, ~] = get_x_multi(N, LV, observations, Theta, deadline, 1);

Y = zeros(N, 3);
Y(:,1) = sum(abs(log(X) - log(observations.y1))<1, 2);
if isfield(observations, 'y2')  
    Y(:,2) = Y(:,1);
    Y(:,3) = sum(abs(log(X) - log(observations.y2))<1, 2);
    Y(:,1) = max(Y(:,1:2)')';
end

keep = (Y(:,3)==10);%(Rej==0);

Theta_out = Theta(keep,:);
X_out = X(keep,:);
n = sum(keep);
end
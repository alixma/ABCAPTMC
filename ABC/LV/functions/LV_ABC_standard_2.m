function [Theta, X, Rej, n, ne] = LV_ABC_standard_2(params, LV, observations, epsilon, cut)
% Standard ABC algorithm for Lotka-Volterra example
nET = observations.ET-1; deadline.T = params.T; N=params.N; %prior = params.prior;
burnin=params.burnin;burning=1; SIGMA = params.SIGMA(:,1);
if~exist('cut', 'var')
    cut=1;
end

% initialise theta and x for each chain
theta = params.theta_in; %exprnd(1./prior);
num_params = max(size(theta));

Theta = zeros(num_params, N);
Rej = zeros(1, N);
X = zeros(nET, N);
Theta(:, 1) = theta;
iter = 1:1000:N;
t=0; ti=1;
deadline.all=tic;
x = get_x_multi(1, LV, observations, theta, deadline);
X(:, 1) = x;

%% LOCAL MOVES %%
while (ti < N)&&(t<=params.T) %local move
    t = toc(deadline.all);
    if burning&&(t>burnin)
        ne = ti;
        burning=0;
    end
    
    if(sum(ti==iter)~=0)
        fprintf('iteration %d, time %f \n', ti, t)
    end
    [theta, x, rej]=LV_local_moves_standard(deadline, theta, x, observations, params, LV, epsilon, SIGMA);
    Theta(:, ti+1) = theta;
    X(:, ti+1) = x;
    Rej(ti+1) = rej;  
    ti = ti+1;
    if(ti == N-100)
        N = N*2;
        iter = fix(linspace(N/2+1, N, 1000));
        fprintf('Extend N to %d \n', N)
    end
end

if(cut)
    Theta = Theta(:,1:ti);
    X = X(:, 1:ti);
    Rej = Rej(:, 1:ti);
end
    n=ti;
end
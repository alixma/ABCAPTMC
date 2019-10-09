function [Theta, X, Rej, n, ne, TM] = LV_ABC_standard(params, LV, observations, epsilon, cut)
% Standard ABC algorithm (1-hit kernel) for Lotka-Volterra example
ET = observations.ET; deadline.T = params.T; 
burning = 1; SIGMA = params.SIGMA(:,1);

if~exist('cut', 'var')
    cut = 1;
end

% initialise theta and x for each chain
theta = params.theta_in; %exprnd(1./params.prior);

Theta = zeros(3, params.N);
Rej = zeros(1, params.N);
X = zeros(ET-1, params.N);
Theta(:, 1) = theta;
iter = 1:100:params.N;
tis = 5:5:params.T+100;
TM = zeros(length(tis), 2); tidx=1;
t=0; ti=1;
deadline.all = tic;
x = get_x_multi(1, LV, observations, theta, deadline);
X(:, 1) = x;

%% LOCAL MOVES %%
while (ti < params.N)&&(t<=params.T) 
    t = toc(deadline.all);
   
    if burning&&(t>params.burnin)
        ne = ti;
        burning = 0;
    end
    
    if(sum(ti==iter)~=0)
        fprintf('iteration %d, time %f \n', ti, t)
    end
    rej = -2;
    while(rej==-2) % try again with new proposal if race took too long
        [theta, x, rej] = LV_local_moves_standard(deadline, theta, x, observations, params, LV, epsilon, SIGMA);
    end
    Theta(:, ti+1) = theta;
    X(:, ti+1) = x;
    Rej(ti+1) = rej;
    ti = ti + 1;
    if(ti == params.N-100)
        params.N = params.N*2;
        iter = fix(linspace(params.N/2+1, params.N, 1000));
        fprintf('Extend params.N to %d \n', params.N)
    end
    if t>tis(tidx) %record sample size
        TM(tidx, :) = [t ti];
        tidx = tidx + 1;
    end
end

if(cut)
    Theta = Theta(:,1:ti);
    X = X(:, 1:ti);
    Rej = Rej(:, 1:ti);
    TM = TM(1:tidx-1,:);
end
n = ti;
end
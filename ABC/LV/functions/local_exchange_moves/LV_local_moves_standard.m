function [theta_out, x_out, rej_out] = LV_local_moves_standard(deadline, theta_in, x_in, observations, params, LV, epsilon, SIGMA)
% Standard local ABC moves for LV model
% Inputs:       deadline: time and deadline parameters
%               theta_in, x_in: input states
%               observations: observations and simulation settings
%               params: various parameters e.g. current state, prior
%               LV: Lotka-Volterra model settings
%               epsilon: tolerance level
%               SIGMA: covariance deviation for simulating from likelihood
%               prior_params: prior parameters
% Outputs:      theta_out, x_out: updated state
%               rej_out: proposal rejected by race(1), prior(2) or not(0)

theta = theta_in; x = x_in;
rej_out = 1; deadline.epsilon = epsilon; SIGMA = diag(SIGMA);

%% PROPOSAL %%
% truncated normal
l = zeros(size(theta));

u = l + observations.run_multi*3 + (1-observations.run_multi)*10; % for single processor example
theta_star = theta + mvrandn(l-theta, u-theta, SIGMA, 1)';

%% PRIOR CHECK %%
A = params.pr(theta_star)*tmvnpdf(theta, theta_star, SIGMA, l, u)/(params.pr(theta)*tmvnpdf(theta_star, theta, SIGMA, l, u));
if(unifrnd(0, 1) < A)
    x_star = zeros(1, observations.ET-1); 
    race = 1; 
else
    race = 0; star=0; rej_out=2;
end

%% RACE %%
rt = tic;
while(race)
    [x, rej] = get_x(LV, observations, theta, deadline);
    [x_star, rej_star] = get_x(LV, observations, theta_star, deadline);

    if(rej+rej_star<0) % deadline was reached
        star = 0; rej_out = -1; %Time out
        break
    end
    
    if(toc(rt)>1000) % 100 % taking too long, start over
        star = 0; rej_out = -2;
        % print('reset\n')
        break
    end 
    
    acc_y = [sum(abs(log(x)-log(observations.y))<=epsilon) sum(abs(log(x_star) - log(observations.y))<=epsilon)];
    acc_te = (acc_y == length(observations.y));
    
    if(sum(acc_te)>=1)
        star = (acc_te(2)>0);
        break
    end
end

%% UPDATE %%
if(star)
    theta = theta_star;
    x = x_star;
    rej_out = 0;
end

theta_out = theta;
x_out = x; 
end


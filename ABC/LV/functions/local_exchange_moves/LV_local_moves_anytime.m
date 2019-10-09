function [theta_out, x_out, resume, resumesim, ti, rej_out] = LV_local_moves_anytime(deadline, observations, params, LV, exchange, epsilon, SIGMA, resume, resumesim)
% Anytime local ABC moves for LV model
% Inputs:       deadline: time and deadline parameters
%               observations: observations and simulation settings
%               params: various parameters e.g. current state, prior
%               LV: Lotka-Volterra model settings
%               epsilon: tolerance level
%               SIGMA: covariance simulating from likelihood
%               prior_params: prior parameters
%               resume: true(1) if a race is to be resumed
%               resumesim: true(1) if a simulation is to be resumed
% Outputs:      theta_out, x_out: updated state
%               resume: true(1) if the race was interrupted
%               resumesim: true(1) if simulation was interrupted
%               ti: time spent in the race
%               rej_out: proposal rejected by race(1), prior(2) or not(0)

theta = params.theta_current(1,:);
x = params.x_current(1,:);
rej_out = 1; deadline.epsilon = epsilon;
SIGMA = diag(SIGMA);

if(resume) % no new proposal if resuming after interruption
    %% RESUMING %%
    theta_star = params.theta_current(2,:);
    theta = params.theta_current(1,:);
    x_star = params.x_current(2,:);
    x = params.x_current(1,:);
    race = 1; % resume race immediately
else    
    %% PROPOSAL %%
    l = zeros(size(theta));
    u = l + 10;
    % use truncated normal for new proposal
    theta_star = theta + mvrandn(l-theta, u-theta, SIGMA, 1)'; 
    x_star = zeros(1, observations.ET-1);
    %% PRIOR CHECK %%    
    A = params.pr(theta_star)*tmvnpdf(theta, theta_star, SIGMA, l, u)/(params.pr(theta)*tmvnpdf(theta_star, theta, SIGMA, l, u));
    if(unifrnd(0, 1) < A)
        race = 1;
    else
        race = 0; star = 0; rej_out = 2;
    end
end

%% RACE %%
while(race)
    deadline.ti = toc(deadline.r);
    if(deadline.ti>exchange.deltat)||((deadline.ti+deadline.t)>params.T) % interrupt the race if either deadline is met
        resume = 1;
        theta_out = [theta; theta_star];
        % fprintf('Stopped at deadline \n')
        x_out = [params.x_current(1,:); x_star];
        break
    end
    
    if(resumesim.resume) % if simulation was interrupted, then resume simulations
        if(resumesim.star) % resume x
            [x, ~, deadline.ti, resumesim] = get_x(LV, observations, theta, deadline, resumesim);
        else % resume x_star
            [x_star, ~, deadline.ti, resumesim] = get_x(LV, observations, theta_star, deadline, resumesim);
        end
        % fprintf('Resume simulations in race \n')
    else % else simulate new x and x_star
        % if interruption was not in simulation then simulate new x and x_star
        resumesim.star = 1; % remember x was interrupted
        [x, ~, deadline.ti, resumesim] = get_x(LV, observations, theta, deadline, resumesim);
        if(~resumesim.resume) % no x_star if simulation on x was interrupted
            resumesim.star=0; % remember x_star was interrupted
            [x_star, ~, deadline.ti, resumesim] = get_x(LV, observations, theta_star, deadline, resumesim);
        end
    end
    resume = resumesim.resume;
    
    if(resumesim.resume) % stopping because simulation was interrupted
        % fprintf('Stopped in race simulations \n')
        theta_out = [theta; theta_star];
        x_out = [params.x_current(1,:); x_star];
        break
    end
    
    if(~resume) % don't compare x and x_star if simulation was interrupted
        acc_y = [sum(abs(log(x)-log(observations.y))<=epsilon) sum(abs(log(x_star) - log(observations.y))<=epsilon)];
        acc_te = acc_y == length(observations.y);
        if(sum(acc_te)>=1)
            % fprintf('Stopped when race ended \n')
            star = acc_te(2)>0;
            break
        end
    end
end

%% UPDATE %%
if(~resume)
    % fprintf('Reached the update step \n')
    if(star)
        theta = theta_star;
        x = x_star;
        rej_out = 0;
    end
    
    theta_out = theta;
    x_out = x;
end
ti = deadline.ti;
end
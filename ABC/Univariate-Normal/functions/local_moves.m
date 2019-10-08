function [theta_out, x_out, resume, ti, rej] = local_moves_adaptive(x_in, theta_in, timepar, params, exchange, observations, epsilon)
% Function to perform local ABC moves
% Inputs:       Tt: {intermediate deadline, 'real' or 'virtual' time}
%               T: global deadline
%               t: to keep track of time until global deadline
%               ti, r: to keep track of time until intermediate deadline
%               y: observations
%               x_in, theta_in: current state
%               sigma: standard deviation for simulating from likelihood
%               epsilon: tolerance level
%               prior_params: prior parameters
%               rho: standard deviation for random walk metropolis proposal
%               pre: include prior check (true (1) by default)
%               resume: true (1) if a race is to be resumed
% Outputs:      theta_out, x_out: updated state
%               resume: true (1) if the race was interrupted
%               ti: time spent in the race 
rej = 1; % proposal lost race
T = timepar.T;
ti = timepar.ti; %t = timepar.t;
resume = timepar.resume;
theta = theta_in(1);
x = x_in(1);

%% PROPOSAL %%
if(resume) %proposal from resuming previous race
    %% RESUMING %%
    theta_star = theta_in(2);
    x_star = x_in(2);    
    race=1;
else %new random walk metropolis proposal
    %% NEW PROPOSAL %%
    theta_star = normrnd(theta, params.rho);
    x_star = 0;
    %% PRIOR CHECK %%
    % using the prior, determine whether theta_star is likely to be rejected
    A = normpdf(theta_star, params.prior(1), params.prior(2))/normpdf(theta, params.prior(1), params.prior(2));
    if(A > unifrnd(0, 1))
        race = 1;
    else
        race = 0; star = 0; 
        rej = 2; % proposal rejected at prior check
    end       
end



%% RACE %%
while(race)
    ti = toc(timepar.r);
        
    if(ti>exchange.deltat)||(toc(timepar.all)>T) % interrupt the race if either deadline is met        
        resume = 1; % resume race after interruption
        theta_out = [theta theta_star];
        x_out = [x_in(1) x_star];
        break
    end

    x = normrnd(theta, observations.sigma);
    x_star = normrnd(theta_star, observations.sigma);
    
    % determine whether x or x_star hit the ball of radius epsilon
    acc_y = [abs(observations.y-x) abs(observations.y-x_star)];    
    acc_te = acc_y<epsilon;    
    acc_tes = sum(acc_te);
    
    if(acc_tes>0)
        %fprintf('race ended, t=%f, ti=%f \n', t, ti);
        star = acc_te(2) > 0;  
        resume = 0;
        break
    end
    
end

%% UPDATE %%
if(~resume)
    if(star)
        rej = 0; %proposal accepted
        theta = theta_star;
        x = x_star;
    end
    
    theta_out = theta;
    x_out = x;
end
end
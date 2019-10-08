function [theta_out, x_out, resume, resumesim, ti, rej_out]=LV_local_moves_anytime(deadline, observations, params, LV, exchange, epsilon, SIGMA, resume, resumesim)
% Anytime local ABC moves for LV model
% Inputs:       Tt: intermediate deadline
%               T: global deadline
%               t: to keep track of time until global deadline
%               ti, r: to keep track of time until intermediate deadline
%               y: observations
%               x_in, theta_in: current state
%               sigma: standard deviation for simulating from likelihood
%               epsilon: tolerance level
%               prior_params: prior parameters
%               resume: true (1) if a race is to be resumed
% Outputs:      theta_out, x_out: updated state
%               resume: true (1) if the race was interrupted
%               ti: time spent in the race

theta = params.theta_current(1,:);
x = params.x_current(1,:);
rej_out = 1; deadline.epsilon = epsilon;
deltat = exchange.deltat; y = observations.y;
pr = params.pr; T = params.T; SIGMA = diag(SIGMA);

if(resume) %no new proposal if resuming after interruption
    %% RESUMING %%
    theta_star = params.theta_current(2,:);
    theta = params.theta_current(1,:);
    x_star = params.x_current(2,:);
    x = params.x_current(1,:);
    race=1; %resume race immediately
else    
%% PROPOSAL %%
    l=zeros(size(theta));
    u=l+10;
    %truncated normal
    theta_star = theta + mvrandn(l-theta, u-theta, SIGMA, 1)'; %new proposal
    x_star = zeros(1, observations.ET-1);
    %% PRIOR CHECK %%    
    A = pr(theta_star)*tmvnpdf(theta, theta_star, SIGMA, l, u)/(pr(theta)*tmvnpdf(theta_star, theta, SIGMA, l, u));
    if(unifrnd(0, 1) < A)
        race = 1;
    else
        race = 0; star=0; rej_out=2;
    end
end

%% RACE %%
while(race)
    
    deadline.ti = toc(deadline.r);
    if(deadline.ti>deltat)||((deadline.ti+deadline.t)>T) %interrupt the race if either deadline is met
        resume = 1;
%         star=0;
        theta_out = [theta; theta_star];
        %fprintf('Stopped at deadline, resume=%f, resumesim=%f \n', resume, resumesim.resume)
        x_out = [params.x_current(1,:); x_star];
        break
    end
    
    if(resumesim.resume) % if simulation was interrupted, then resume simulations
        if(resumesim.star) %resume x
            [x, ~, deadline.ti, resumesim] = get_x(LV, observations, theta, deadline, resumesim);
        else %resume x_star
            [x_star,~,deadline.ti, resumesim] = get_x(LV, observations, theta_star, deadline, resumesim);
        end
        %fprintf('Resume simulations in race, resumesim=%d resume=%d \n', resumesim.resume, resume)
    %end % ADDED 08/02/2018
    else % ADDED 11/05/18 else simulate new x and x_star
        % if interruption was not in simulation then simulate new x and x_star
        resumesim.star=1;
        [x, ~, deadline.ti, resumesim] = get_x(LV, observations, theta, deadline, resumesim);
        if(~resumesim.resume) %EDITED 08/02/2018 %no x_star if simulation on x was interrupted
            resumesim.star=0; %remember x_star was interrupted
            [x_star,~,deadline.ti, resumesim] = get_x(LV, observations, theta_star, deadline, resumesim);
        end
    end
    resume=resumesim.resume;
    
    if(resumesim.resume) %stopping because simulation was interrupted
        %resume=1;
        %star=0;
        %fprintf('Stopped in race simulations, resume=%f, resumesim=%f \n', resume, resumesim.resume)
        theta_out = [theta; theta_star];
        x_out = [params.x_current(1,:); x_star];
        break
    end
    
    if(~resume) %EDITED 08/02/2018 don't compare x and x_star if simulation was interrupted
        acc_y = [sum(abs(log(x)-log(y))<=epsilon) sum(abs(log(x_star) - log(y))<=epsilon)];
        acc_te = acc_y == length(y);
        if(sum(acc_te)>=1)
            %fprintf('Stopped in race, resume=%f, resumesim=%f \n', resume, resumesim.resume)
            star=acc_te(2)>0;
            break
        end
    end
end

%% UPDATE %%
%resume
%resumesim.resume
if(~resume)
    %fprintf('Reached the last resume, resume=%f, resumesim=%f \n', resume, resumesim.resume)
    if(star)
        theta = theta_star;
        % params.x_current = x_star;  EDIT 09/01/2019
        x = x_star;
        %EDIT 09/01/2019
        rej_out = 0;
    end
    
    theta_out = theta;
    %x_out = params.x_current(1,:); %EDIT 09/01/2019: if proposal rejected, keep x that won the
    %race rather than initial x.    
    x_out = x;
    %EDIT 09/01/2019
end


ti = deadline.ti;%toc(deadline.r);

end
function [theta_out, x_out, rej_out]=LV_local_moves_standard(deadline, theta_in, x_in, observations, params, LV, epsilon, SIGMA)
% Standard (not anytime) local ABC moves for LV model
deadline.epsilon = epsilon;
SIGMA = diag(SIGMA); 
theta = theta_in;
x = x_in;
rej_out = 1;

%% PROPOSAL %%
%truncated normal
l=zeros(size(theta));
u=l+3;
theta_star = theta + mvrandn(l-theta, u-theta, SIGMA, 1)';

%% PRIOR CHECK %%
A = params.pr(theta_star)*tmvnpdf(theta, theta_star, SIGMA, l, u)/(params.pr(theta)*tmvnpdf(theta_star, theta, SIGMA, l, u));
if(unifrnd(0, 1) < A)
    x_star = zeros(1, observations.ET-1); %get_x(LV, observations, theta_star, deadline);
    race = 1; %resume = 1;
else
    race = 0; star=0; rej_out=2;
end


%% RACE %%
rt = tic; % AEDDD 20/09/19 
while(race)
    [x, rej] = get_x(LV, observations, theta, deadline);
    [x_star, rej_star] = get_x(LV, observations, theta_star, deadline);

    if(rej+rej_star<0) %deadline was reached
        star=0; rej_out=-1; %Time out
        break
    end
    
    if(toc(rt)>5) %taking too long, start over
        star=0; rej_out=-2; %Time out
        %print('reset\n')
        break
    end 
    
    acc_y = [sum(abs(log(x)-log(observations.y))<=epsilon) sum(abs(log(x_star) - log(observations.y))<=epsilon)];
    acc_te = (acc_y == length(observations.y));
    
%     if isfield(observations, 'y2')
%         acc_y2 = [sum(abs(log(x)-log(observations.y2))<=epsilon) sum(abs(log(x_star) - log(observations.y2))<=epsilon)];
%         acc_te2 = (acc_y2 == length(observations.y2));
%     else
%         acc_te2 = zeros(1, 2);
%     end
    if(sum(acc_te)>=1)%||(sum(acc_te2)>=1)
        star=(acc_te(2)>0);%||(acc_te2(2)>0);
        break
    end
end

%% UPDATE %%
if(star)
    theta = theta_star;
    %x_in = x_star; EDIT 09/01/2019
    x = x_star; %
    %/EDIT 09/01/2019
    rej_out = 0;
end

theta_out = theta;
%x_out = x_in; %EDIT 09/01/2019: if proposal rejected, keep x that won the
%race rather than initial x.
x_out = x; 
%/EDIT 09/01/2019

end


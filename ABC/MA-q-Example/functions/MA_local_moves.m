function [theta_out, x_out, resume, ti, rej_out] = MA_local_moves(r, all, ti, theta_in, x_in, params, rho, simulations, epsilon, resume)

rej_out = 1;
%% PROPOSAL %%
if(resume) %proposal from resuming previous race
    theta_star = theta_in(:,2);
    theta = theta_in(:,1);
    x_star = x_in(:,2);
    x = x_in(:,1);
    race = 1; %resume race
else % propose new thetas (random walk proposal)
    theta = theta_in(:,1);
    x = x_in(:,1);
    theta_star = normrnd(theta, rho);
    x_star = 10+zeros(simulations.q, 1);%MA_get_stats(theta_star, simulations.q, nsamples);
    
%% PRIOR CHECK
% reject if outside the triangle    
    if(theta_star(2)<=-1)||(theta_star(2)>=1)||(theta_star(1)<=-2)||(theta_star(1)>=2)||(sum(theta_star)<=-1)||(theta_star(1)-theta_star(2)>=1)
        star=0;
        race=0; %do not race if outside triangle
        rej_out = 2;
    else
        race=1; %race if inside triangle
    end
end
%E = epsilon(min(max(1,floor(toc(all))), params.T));


%% RACE %%
while(race)
    ti = toc(r);
    
    % interrupt race if either deadline is met
    if((ti>params.deltat) ||(toc(all)>params.T))
        resume=1;
        star=0;
        theta_out = [theta, theta_star];
        x_out = [x_in(:,1), x_star];
        break
    end

    
    x = MA_get_stats(theta, simulations);
    x_star = MA_get_stats(theta_star, simulations);
    resume=0;

    
    acc_y = [norm(simulations.y'-x), norm(simulations.y'-x_star)];
    acc_te = acc_y < epsilon;
    
    if(sum(acc_te)>0)
        star = acc_te(2)>0;
        %race=0;
        break
    end
end
%   mess='reject';
if(~resume)
    if(star)
%           mess = 'accept';
        theta = theta_star;
        x = x_star;
        rej_out = 0;
    end
    theta_out = theta;
    x_out = x; %x_in(:,1);
%       fprintf('%s proposal, distance=%f, epsilon=%f \n\n', mess, norm(simulations.y'-x_out), E);
end

end
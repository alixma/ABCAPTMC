function [xmat, rej, ti, resumesim] = gillespiedy_anytime(LV, observations, th, deadline, resumesim)
% Discretised Gillespie algorithm for a SPN
% WITH DEADLINE controlled by TIMER

% initialise system at t = 0 (or resume)
n = fix(observations.ET/observations.dt);
S=(LV.Post-LV.Pre)';
[u, v] = size(S);
rej= 0; repeat = 1;

if(resumesim.resume)
    target = resumesim.target; xmat = resumesim.xmat;
    i = resumesim.i; evolve = resumesim.evolve;
    x = resumesim.x; tt = resumesim.tt;
else    
    tt = 0; target = 0;
    xmat = zeros(n, u); i = 1;
    evolve = 1; x = LV.M;
end

while(repeat)
    % deadline checkpoint
    ti = toc(deadline.r);
    if(ti>deadline.deltat)||((ti+deadline.t)>deadline.T)
        % fprintf('Time is up \n');
        resumesim.xmat=xmat; resumesim.target = target;
        resumesim.resume=1; resumesim.i=i; % simulation will be resumed
        resumesim.x=x; resumesim.tt= tt; resumesim.evolve = evolve;
        break
    end    
    % compute hazard function
    h = LV.h(x, th); h0 = sum(h);  
    if(h0<1e-10)
        % fprintf('Extinct! \n')
        evolve = 0;
        tt = 1e99;
    else        
        % simulate time to next event based on combined reaction hazard
        tt = tt + exprnd(1/h0);
    end
    
    % when intermediate (discrete) target is reached
    while(tt>=target)
        if deadline.early
            epsilon_check = ((target>0)&&(abs(log(x(1)) - log(observations.y(target)))>deadline.epsilon));
        end
        if deadline.early&&(h0>1e3)&&epsilon_check
            % if prey will not hit the ball or time is up
            % fprintf('prey will not hit ball, target=%d, x=%d, y=%d, i=%d \n', target, x(1), y(target), i);
            xmat(i,:) = x;
            rej=1;
            resumesim = []; resumesim.resume = 0; % simulation will not be resumed
            repeat=0; % abort
            break
        else
            %record x
            xmat(i,:) = x;
            i = i+1;
            target = target + observations.dt;
            if(i>n)
                repeat = 0; 
                resumesim = []; resumesim.resume = 0; %simulation will not be resumed
                break
            end
        end
    end
    
    % population evolves
    if(repeat&&evolve)
        % simulate reaction index if not extinct:
        %   j=1 prey birth
        %   j=2 prey consumption
        %   j=3 predator death
        j = randsample(v, 1, true, h);
        % update x according to reaction
        x = x + S(:,j);
        if((x(2)==0)&&(target<observations.ET-3*observations.dt)) % if predator goes extinct
            % fprintf('predator Extinct! \n')
            xmat(:,1) = 1e12; % ensure xmat is rejected
            rej=1+sum(find(x==0));
            repeat=0; resumesim=[]; resumesim.resume=0; % simulation will not be resumed
        end
    end
end

end
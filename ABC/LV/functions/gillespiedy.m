function [xmat, rej] = gillespiedy(LV, observations, th, deadline)%, y, epsilon)
% Discretised Gillespie algorithm for a SPN (TIMER DEADLINE)



% initialise system at t=0
n = fix(observations.ET/observations.dt);
S=(LV.Post-LV.Pre)';
rej=0;
[u, v] = size(S);
tt=0;
target=0;
xmat = zeros(n, u);
i=1;
x=LV.M;
%tvec=target:dt:(ET-dt);







repeat=1;
evolve=1;

while(repeat)
    % deadline checkpoint
    if(toc(deadline.all)>deadline.T)
        %fprintf('Time is up \n');
        xmat(i,:) = x;
        rej=-1;
        break
    end
    
    
    
    
    
    %compute hazard function
    h=LV.h(x, th);
    h0=sum(h);
    if(h0<1e-10)
        %fprintf('Extinct! \n')
        evolve=0;
        tt = 1e99;
    else
        % simulate time to next event based on combined reaction hazard
        tt = tt + exprnd(1/h0);
    end
    
    % when intermediate (discrete) target is reached
    while(tt>=target)
        if isfield(observations, 'y2')&&deadline.early
            epsilon_check = ((target>0)&&(abs(log(x(1)) - log(observations.y(target)))>deadline.epsilon))&&(abs(log(x(1)) - log(observations.y2(target)))>deadline.epsilon);
        elseif deadline.early
            epsilon_check = ((target>0)&&(abs(log(x(1)) - log(observations.y(target)))>deadline.epsilon));
        end
        if deadline.early&&(h0>1e3)&&epsilon_check
            %if prey will not hit the ball
            %fprintf('prey will not hit ball, target=%d, x=%d, y=%d, i=%d \n', target, x(1), y(target), i);
            %xmat(:,1) = 1e9; % ensure xmat is rejected
            xmat(i,:) = x;
            rej=1;
            repeat=0; %abort
            break
        else
            %record x
            xmat(i,:) = x;
            i = i+1;
            target = target+observations.dt;
            if(i>n)
                repeat=0;
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
        if((x(2)==0)&&(target<observations.ET-3*observations.dt))%||(x(1)==0)%% %if predator (EDIT:or prey, commented out) goes extinct
            %fprintf('predator Extinct! \n')
            xmat(:,1) = 1e12; % ensure xmat is rejected
            rej=1+sum(find(x==0));
            repeat=0; %abort
        end
    end
end

end
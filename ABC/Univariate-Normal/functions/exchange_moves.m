function [theta_out, x_out, n_in, rej_out] = exchange_moves(theta_in, x_in, n_in, swap, epsE, tp, observations, simple)
rej_out = 5; % indicates 'race' not completed
theta_out = theta_in;
x_out = x_in;

if(simple)
    % simple
    %% MULTIPLE SIMPLE SWAPS
    small = swap(:,1); E_small = epsE(small)';
    big = swap(:,2); %E_big = epsE(big)';
    acc_y = abs(x_in(big)-observations.y);
    mkswaps = acc_y<E_small;
    
    % perform exchange moves where accepted
    theta_out(big(mkswaps)) = theta_in(small(mkswaps)); 
    theta_out(small(mkswaps)) = theta_in(big(mkswaps)); 
    x_out(big(mkswaps)) = x_in(small(mkswaps)); 
    x_out(small(mkswaps)) = x_in(big(mkswaps)); 
    
    n_in(swap(:)) = n_in(swap(:))+1;  
    rej_out = 3 + mkswaps; % 4 for accepted 3 for rejected swaps
else
    % new 'race' before performing exchange moves
    %% ONE NOT SO SIMPLE SWAP
    sample = 1;
    E = epsE(swap);
    while(sample)
        x_star = normrnd(theta_in(swap(2)), observations.sigma);
        acc_y = abs(x_star-observations.y);
        tp.t=toc(tp.all); % check time
        if(tp.t>tp.T)  % stop if global deadline occurs
            mkswap = 0;
            rej_out=5;
            break
        end             
        % stop if larger ball is hit
        if(acc_y<E(2)||tp.t>tp.T)
            mkswap = tp.t<=tp.T;
            break;
        end
    end
    % attempt to make swap if hard deadline hasn't been reached
    if(mkswap)
        rej_out = 3; % indicates 'race' completed but swap failed
        x_in(swap(2)) = x_star;
        %accept if smaller ball is hit
        if(acc_y<E(1))
            rej_out = 4; % indicates swap succeeded           
            % perform the accepted swaps
            theta_out(swap(2)) = theta_in(swap(1)); x_out(swap(2)) = x_in(swap(1));
            theta_out(swap(1)) = theta_in(swap(2)); x_out(swap(1)) = x_in(swap(2));
        end        
        n_in(swap) = n_in(swap)+1;        
    end
end
end
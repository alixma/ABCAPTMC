function [theta_out, x_out, n_out, rej_out] = MA_exchange_moves(theta_in, x_in, n_in, simulations, epsilon, swap, sz, fast)
% Exchange moves for MA example
            
if ~exist('fast', 'var')
    fast=1;
end

theta_out = theta_in; x_out = x_in; n_out = n_in;

if fast % fast exchange moves
    if ~exist('sz', 'var')
        sz = size(swap, 1);
    end
    rej_out = 3+zeros(sz, 1); %3 if rejected swap
    
    % check whether smaller ball is also hit
    acc_y = vecnorm(x_in(:, swap(:,2))-simulations.y');
    mkswap = acc_y < epsilon((swap(:,1)));
    
    % perform accepted swaps
    rej_out(mkswap) = 4; %accepted swaps
    theta_out(:, swap(mkswap, 2)) = theta_in(:, swap(mkswap, 1)); x_out(:, swap(mkswap, 2)) = x_in(:, swap(mkswap, 1));
    theta_out(:, swap(mkswap, 1)) = theta_in(:, swap(mkswap, 2)); x_out(:, swap(mkswap, 1)) = x_in(:, swap(mkswap, 2));
else %  slow exchange moves
    rej_out = 3; sample = 1;
    % sample until larger ball is hit
    while(sample)

        x_star = MA_get_stats(theta_in(:,swap(2)), simulations);
        acc_y = norm(x_star-simulations.y');
        t=toc(all);
        if(t>params1.T)
            mkswap=0;
            break
        end
        E = epsilon(swap);

        if(acc_y<E(2))||(t>params1.T)
            mkswap=t<=params1.T;
            sample=0;
        end
    end
    % attempt swap (unless global deadline occurred)
    if(mkswap)
        %             fprintf('Worker %d: Attempting swap t= %f \n', c, t)
        x_in(:,swap(2)) = x_star;
        %accept or reject swap
        if(acc_y<E(1))
            rej_out = 4; %accepted swap
            % perform the accepted swaps
            theta_out(:, swap(2)) = theta_in(:, swap(1)); x_out(:, swap(2)) = x_in(:, swap(1));
            theta_out(:, swap(1)) = theta_in(:, swap(2)); x_out(:, swap(1)) = x_in(:, swap(2));
        end

    end
end
n_out(swap) = n_in(swap)+1;   
end
    
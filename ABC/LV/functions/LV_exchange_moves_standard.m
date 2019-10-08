function [theta_out, x_out, rej_out, n_out, swap]=LV_exchange_moves_standard(deadline, theta_in, x_in, n_in, observations, LV, exchange, epsilon, swap)
% Standard (not anytime) exchange ABC moves for LV model

sample=1; rej_out=3; nET = observations.ET-1;
theta_out = theta_in; x_out = x_in; n_out = n_in;
if ~exist('swap', 'var')
    swap = sort(exchange.pair(randsample(exchange.ipairs, 1), :));
end
deadline.epsilon = epsilon(swap(2)); % for early interruption of simulation



% sample until larger ball is hit
% while(sample)
%     [x_swap, rej_swap] = get_x(LV, observations, theta_in(swap(2),:), deadline);
%     acc_y = abs(log(x_swap)-log(observations.y));
% 
%     if(sum(acc_y<epsilon(swap(2)))==nET)||(rej_swap==-1)
%         mkswap=(sum(acc_y<epsilon(swap(1)))==nET);
%         break
%     end
% end
% fast exchanges
acc_y = abs(log(x_in(swap(2),:))-log(observations.y));
mkswap=(sum(acc_y<epsilon(swap(1)))==nET);

% make swap (if smaller ball was hit)
if(mkswap)
    rej_out = 4; %accepted swap
    %x1(swap_w(2), :) = x_swap;
    theta_out(swap(2),:) = theta_in(swap(1),:); x_out(swap(2),:) = x_in(swap(1),:);
    theta_out(swap(1),:) = theta_in(swap(2),:); x_out(swap(1),:) = x_in(swap(2),:);
end
n_out(swap) = n_in(swap)+1;
end


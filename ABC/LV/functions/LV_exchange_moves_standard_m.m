function [theta_out, x_out, rej_out, n_out, swap]=LV_exchange_moves_standard_m(theta_in, x_in, n_in, observations, epsilon, sz, swap)
% Fast exchange ABC moves for LV model, multiple at once

if ~exist('sz', 'var')
        sz = size(swap, 1);
end
rej_out=3+zeros(sz, 1); nET = observations.ET-1;
theta_out = theta_in; x_out = x_in; n_out = n_in;

acc_y = abs(log(x_in(swap(:, 2),:))-log(observations.y));
mkswap=(sum(acc_y<repmat(epsilon(swap(:, 1))', 1, nET), 2)==nET);

% make swaps (if smaller balls was hit)
rej_out(mkswap) = 4; %accepted swaps
theta_out(swap(mkswap, 2), :) = theta_in(swap(mkswap, 1),:); x_out(swap(mkswap, 2),:) = x_in(swap(mkswap, 1),:);
theta_out(swap(mkswap, 1), :) = theta_in(swap(mkswap, 2),:); x_out(swap(mkswap, 1),:) = x_in(swap(mkswap, 2),:);

n_out(swap) = n_in(swap)+1;
end


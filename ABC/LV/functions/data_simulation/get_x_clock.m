function [x_out, rej_out, res_out] = get_x_clock(LV, observations, th, deadline, res_in)%, show)
% Function to obtain simulated values of X_1 at times 1,...10
deadline.early=1;
if~exist('res_in', 'var')
    [x, rej_out] = gillespiedy_clock(LV, observations, th, deadline);
    x_out = x(2:end, 1)';
else
    [x, rej_out, res_out] = gillespiedy_anytime_clock(LV, observations, th, deadline, res_in);
    x_out = x(2:end, 1)';
end
end
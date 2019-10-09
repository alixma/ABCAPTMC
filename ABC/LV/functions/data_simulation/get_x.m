function [x_out, rej_out, ti, res_out]=get_x(LV, observations, th, deadline, res_in)
% Function to obtain simulated values of X_1 at times 1,...10

deadline.early = 0;
if~exist('res_in', 'var')
    [x, rej_out] = gillespiedy(LV, observations, th, deadline);
    x_out = x(2:end, 1)';
else
    [x, rej_out, ti, res_out] = gillespiedy_anytime(LV, observations, th, deadline, res_in);
    x_out = x(2:end, 1)';
end

end
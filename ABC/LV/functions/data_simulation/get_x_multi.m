function [x_out, rej_out, ti, res_out]=get_x_multi(K, LV, observations, th, deadline, res_in)
% Function to obtain simulated values of X_1 at times 1,...10
deadline.early=0;

if~exist('res_in', 'var')
    x_out = zeros(K, observations.ET-1);
    rej_out = zeros(K,1);
    for k=1:K
        [x, rej] = gillespiedy(LV, observations, th(k,:), deadline);
        x_out(k,:) = x(2:end, 1);
        rej_out(k) = rej;
    end
else
    resume=res_in;
    x_out = zeros(K, observations.ET-1);
    rej_out = zeros(K,1);
    for k=1:K
        [x, rej, ti, resume] = gillespiedy_anytime(LV, observations, th(k,:), deadline, resume);
        x_out(k,:) = x(2:end, 1);
        rej_out(k) = rej;
        res_out=resume;       
    end
end

end
function [x_out, rej_out, ti, res_out]=get_x_multi(K, LV, observations, th, deadline, res_in)%, show)
% Function to obtain simulated values of X_1 at times 1,...10
if(~exist('show', 'var'))
    show=0;
end
% if(~exist('res_in', 'var'))
%     res_in=0;
% end
deadline.early=0;

if~exist('res_in', 'var')
    x_out = zeros(K, observations.ET-1);
    rej_out = zeros(K,1);
    %ti=NaN; res_out = NaN;
    for k=1:K
%         if(show)&&(mod(k, K/20)==0)
%             fprintf('iteration %d of %d \n', k, K)
%         end
        [x, rej] = gillespiedy(LV, observations, th(k,:), deadline);%, y, epsilon(k));
        x_out(k,:) = x(2:end, 1);
        rej_out(k) = rej;
    end
else
    resume=res_in;
    x_out = zeros(K, observations.ET-1);
    rej_out = zeros(K,1);
    for k=1:K
%         if(show)&&(mod(k, K/20)==0)
%             fprintf('iteration %d of %d \n', k, K)
%         end
        %[x, rej] = gillespiedy(LV, ET, dt, th(k,:));
        [x, rej, ti, resume] = gillespiedy_anytime(LV, observations, th(k,:), deadline, resume);
        x_out(k,:) = x(2:end, 1);
        rej_out(k) = rej;
        res_out=resume;       
    end
    %ti = toc(deadline.r);
end

end
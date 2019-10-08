function [x_out, rej_out, ti, res_out]=get_x(LV, observations, th, deadline, res_in)%, show)
% Function to obtain simulated values of X_1 at times 1,...10
%ET = observations.ET; dt = observations.dt; %y=observations.y;
% if(~exist('show', 'var'))
%     show=0;
% end
% if(~exist('res_in', 'var'))
%     res_in=0;
% end
deadline.early=0;
if~exist('res_in', 'var')%)(nargin<5)
    [x, rej_out] = gillespiedy(LV, observations, th, deadline);%(k,:));%, y, epsilon(k));
    x_out = x(2:end, 1)';
else
%     resume=res_in;    
    [x, rej_out, ti, res_out] = gillespiedy_anytime(LV, observations, th, deadline, res_in);
    x_out = x(2:end, 1)';
%     res_out=resume;    
end

end
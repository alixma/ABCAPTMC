function [ESS_OUT, IAT_OUT] = ess_iat_maxlags(S, maxlag, imse)
if ~exist('imse', 'var')
    imse=1;
end
if ~exist('maxlag', 'var')
    maxlag=1000;
end
ntimes = max(size(S));
ESS_OUT = zeros(ntimes, 3); IAT_OUT = zeros(ntimes, 3);

for w=1:ntimes
    [ESS_OUT(w,:), IAT_OUT(w,:)] = ess_iat_multi({S{w}}, min(maxlag, length(S{w})-1), imse, 0);
end
end
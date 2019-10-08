function [ESS_OUT, IAT_OUT] = autocorr_new_multi(S, maxc)
if ~exist('maxc', 'var')
    maxc=5;
end
ESS_OUT = zeros(1, 3); IAT_OUT = zeros(1, 3);

for w=1:3
    [ESS_OUT(w), IAT_OUT(w)] = autocorr_new(S, w, maxc);
end
end
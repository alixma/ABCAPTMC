function [ess_out, iat_out]=ess_iat_multi(S, maxlag, imse, combine)
% Computing ESS and IAT for W runs
if~exist('combine', 'var')
    combine=0;
end

W = max(size(S)); nparams = size(S{1},2);
ess = zeros(W, nparams); iat = zeros(W, nparams);
for w=1:W
    S_w = S{w};
    for n=1:nparams
        [ess(w, n), iat(w, n)] = ESS_IAT(S_w(:,n), maxlag, imse);
        
    end
end

if(combine)
    ess_out = zeros(1, nparams); iat_out = zeros(1, nparams);
    for n=1:nparams
        ess_out(n) = sum(ess(:,n)); iat_out(n) = mean(iat(:,n));
    end
else
    ess_out = ess;
    iat_out = iat;
end
end


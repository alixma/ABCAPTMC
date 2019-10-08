function [ess_out, iat_out] = ess_iat_allchains(Theta, ne, n, imse, maxlag)
ntimes = size(Theta, 2);
if ~exist('imse', 'var')
    imse=1;
end
if ~exist('maxlag', 'var')
    maxlag=1000;
end
nchains = size(Theta{1}, 1);
S = cell(1, nchains);
for k= 1:nchains     
    S1 = cell(1, ntimes); 
    for idx=1:ntimes
        Theta1 = Theta{idx};
        T1  = Theta1(k, :, ne(k, idx):n(k, idx));
        S1{idx} = reshape(T1, size(T1, 2), size(T1, 3))';
    end
    S{k} = S1;
end
ess_out = zeros(nchains, ntimes, 3); iat_out = zeros(nchains, ntimes, 3);
for k=1:nchains
    [ess_out(k,:,:), iat_out(k,:,:)] = ess_iat_maxlags(S{k}, maxlag, imse);%, 1);
end

end
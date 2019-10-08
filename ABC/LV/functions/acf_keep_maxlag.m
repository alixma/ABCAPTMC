function [R_out, n_algs, Cc]=acf_keep_maxlag(R_in, bs, maxlag)
n_algs = max(size(bs));
Cc = 1e9;
for na=1:n_algs
    Cc = min(Cc, sum(bs{na}>maxlag));
end
R_out = cell(1, n_algs);
s = zeros(n_algs, 1);
for na=1:n_algs
    Cc = sum(bs{na}>maxlag);
    R = cell(1, Cc);
    ntimes = max(size(R_in{na}));
    for ni=1:ntimes
        
        if bs{na}(ni)>maxlag
            s(na)=s(na)+1;
            R{s(na)} = R_in{na}{ni};
        end
        
    end
    R_out{na} = R;
end
end

function save_all(dirname, name, output, n, TM, out, rates)
% save everything
ntimes = max(size(output));

for idx=1:ntimes    
    save_one(idx, dirname, name, output, n, TM, out, rates)
end
end
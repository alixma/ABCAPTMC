function save_all(dirname, name, output, n, TM, out, rates, num)
% save everything
ntimes = max(size(output));
if ~exist('num','var')
    num = 0;
end

for idx=1:ntimes    
    save_one(idx, dirname, name, output, n, TM, out, rates, num)
end
end
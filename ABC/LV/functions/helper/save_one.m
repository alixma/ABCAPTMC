function save_one(idx, dirname, name, output, n, TM, out, rates, num)

if ~exist('num','var')
    num = 0;
end
%output
file = sprintf('results/LV/%s/W%d_%s_output.csv', dirname, idx+num, name);
dlmwrite(file, output{idx});

%sample sizes
file = sprintf('results/LV/%s/W%d_%s_n.csv', dirname, idx+num, name);
dlmwrite(file, n(:,idx));

%timelines
file = sprintf('results/LV/%s/W%d_%s_TM.csv', dirname, idx+num, name);
dlmwrite(file, TM{idx});

%composition of chain
file = sprintf('results/LV/%s/W%d_%s_stats.csv', dirname, idx+num, name);
dlmwrite(file, out(:,:,idx));

%rates
file = sprintf('results/LV/%s/W%d_%s_rates.csv', dirname, idx+num, name);
dlmwrite(file, rates(:,:,idx));

end
function save_one(idx, dirname, name, output, n, TM, out, rates)
%output
file = sprintf('results/LV/%s/W%d_%s_output.csv', dirname, idx, name);
dlmwrite(file, output{idx});

%sample sizes
file = sprintf('results/LV/%s/W%d_%s_n.csv', dirname, idx, name);
dlmwrite(file, n(:,idx));

%timelines
file = sprintf('results/LV/%s/W%d_%s_TM.csv', dirname, idx, name);
dlmwrite(file, TM{idx});

%composition of chain
file = sprintf('results/LV/%s/W%d_%s_stats.csv', dirname, idx, name);
dlmwrite(file, out(:,:,idx));

%rates
file = sprintf('results/LV/%s/W%d_%s_rates.csv', dirname, idx, name);
dlmwrite(file, rates(:,:,idx));

end
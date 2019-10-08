function [output, n, TM, stats, rates] = read_one(dirname, name, ntimes)
% save everything

output = cell(1, ntimes); TM = output;

n1 = dlmread(sprintf('results/ABC/LV/%s/W%d_%s_n.csv', dirname, 1, name));
stats1 = dlmread(sprintf('results/ABC/LV/%s/W%d_%s_stats.csv', dirname, 1, name));
rates1 = dlmread(sprintf('results/ABC/LV/%s/W%d_%s_rates.csv', dirname, 1, name));

n = zeros([max(size(n1)) ntimes]);
stats = zeros([size(stats1) ntimes]);
rates = zeros([size(rates1) ntimes]);
for idx=1:ntimes
    
    %output
    output{idx} = dlmread(sprintf('results/ABC/LV/%s/W%d_%s_output.csv', dirname, idx, name));
    %timelines
    TM{idx} = dlmread(sprintf('results/ABC/LV/%s/W%d_%s_TM.csv', dirname, idx, name));
    
    if idx==1
        n(:, idx) = n1;
        stats(:, :, idx) = stats1;
        rates(:, :, idx) = rates1;
    else
        %sample sizes
        n(:, idx) = dlmread(sprintf('results/ABC/LV/%s/W%d_%s_n.csv', dirname, idx, name));
        
        %composition of chain
        stats(:, :, idx) = dlmread(sprintf('results/ABC/LV/%s/W%d_%s_stats.csv', dirname, idx, name));
        
        %rates
        rates(:, :, idx) = dlmread(sprintf('results/ABC/LV/%s/W%d_%s_rates.csv', dirname, idx, name));
    end
end

end
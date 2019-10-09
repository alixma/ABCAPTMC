function plot_multi_densities(R_in, S_rej, den, names, col)
% Plot marginal densities of three algorithm outputs for comparison
% Algorithms:   1) Standard AMC
%               2) Single processor PTMC
%               3) Single processor APTMC

if ~exist('den', 'var')
    den = [2^8 2^16 2^8
        2^6 2^8 2^6
        2^8 2^8 2^8
        2^6 2^6 2^6];
end

if ~exist('names', 'var')
    names = {'ABC-PTMC-1', 'ABC-APTMC-1', 'Standard ABC'};
end

if ~exist('col', 'var')
    col = {[1 17 181] ./ 255, [251 111 66] ./ 255, [12 195 82] ./ 255};
end
nalgs = size(R_in, 2);

figure;
for i=1:3
    subplot(1, 3, i); hold on; grid on;
    [~,density,xmesh,~] = kde(S_rej(:, i), den(1,1), 0, 4); plot(xmesh, density, 'LineWidth', 1, 'Color',rgb('Black'), 'DisplayName', 'Reference');
    for j=1:nalgs
        [~,density1,xmesh1,~] = kde(R_in{j}(:, i), den(2,1), 0, 4); plot(xmesh1, density1, 'Color',col{j}, 'DisplayName', names{j});
    end
    xlabel(sprintf('\\theta_%d', i));
end

lgd = legend('show');
lgd.FontSize = 7;
lgd.Location = 'north east';
end
function plot_multi_chains(R, true_theta, titl, col)
% Plot MCMC chains output from an ABC-MCMC algorithm
figure; 
W = size(R, 2);
if ~exist('col', 'var')
    col = [rgb('DarkBlue'); rgb('Tomato'); rgb('DarkGreen')];
end
th = length(true_theta);

for t=1:th
    subplot(th, 1, t); hold on;
    for w=1:W
        R_w = R{w}; plot(R_w(:, t));
    end
    title(titl); ylabel(sprintf('\\theta_%d', t));
    hline = refline(0, true_theta(t)); hline.Color = col(t,:);
end

end
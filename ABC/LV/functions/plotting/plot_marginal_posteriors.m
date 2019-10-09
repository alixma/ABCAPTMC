function plot_marginal_posteriors(R, S_rej, bins, titl)
% Plot histograms of the marginal posteriors output from an ABC-MCMC algorithm

figure; hold on;

subplot(1, 3, 1); h=histogram(R(:, 1), 'Normalization', 'pdf');
[~,density,xmesh,~]=kde(S_rej(:, 1), 2^8, 0, 4); hold on;
plot(xmesh, density);
xlabel('\theta_1'); h.NumBins=bins; xlim([-.1 4]);

subplot(1, 3, 2); h=histogram(R(:, 2), 'Normalization', 'pdf');
[~,density,xmesh,~]=kde(S_rej(:, 2), 2^16, 0, 0.04); hold on;
plot(xmesh, density); title(titl);
h.FaceColor = rgb('Orange'); xlabel('\theta_2'); h.NumBins=bins; xlim([-.001 0.04]);

subplot(1, 3, 3); h=histogram(R(:, 3), 'Normalization', 'pdf');
[~,density,xmesh,~]=kde(S_rej(:, 3), 2^8, 0, 4); hold on;
plot(xmesh, density);
h.FaceColor=rgb('Green'); xlabel('\theta_3'); h.NumBins=bins; xlim([-.1 4]);
end
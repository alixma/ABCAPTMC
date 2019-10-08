function plot_chains(R, true_theta, titl)
% Plot MCMC chains output from an ABC-MCMC algorithm
figure; hold on
subplot(3, 1, 1); plot(R(:, 1)); ylabel('\theta_1'); title(titl);
hline=refline(0, true_theta(1)); hline.Color = rgb('DarkBlue');
subplot(3, 1, 2); plot(R(:, 2), 'Color', rgb('Orange')); ylabel('\theta_2');
hline=refline(0, true_theta(2)); hline.Color = rgb('Tomato');
subplot(3, 1, 3); plot(R(:, 3), 'Color', 'g'); ylabel('\theta_3');
hline=refline(0, true_theta(3)); hline.Color = rgb('DarkGreen');
end
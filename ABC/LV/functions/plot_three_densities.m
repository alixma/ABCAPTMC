function plot_three_densities(R1, R2, R3, S_rej, den, names)
% Plot sample acf of three algorithm outputs for comparison
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
    names = {'ABC-PTMC-1', 'ABC-APTMC-1', 'ABC'};
end



colour_darkblue = [1 17 181] ./ 255;
colour_peach = [251 111 66] ./ 255;
colour_green = [12 195 82] ./ 255;

figure;
subplot(1, 3, 1); hold on; grid on;
[~,density,xmesh,~]=kde(S_rej(:, 1), den(1,1), 0, 4); plot(xmesh, density, 'LineWidth', 1, 'Color',rgb('Black'), 'DisplayName', 'Reference');%[0 0.4470 0.7410]);
[~,density1,xmesh1,~]=kde(R1(:, 1), den(2,1), 0, 4); plot(xmesh1, density1, 'Color',colour_darkblue, 'DisplayName', names{1});
[~,density2,xmesh2,~]=kde(R2(:, 1), den(3,1), 0, 4); plot(xmesh2, density2, 'Color',colour_peach, 'DisplayName', names{2});
[~,density3,xmesh3,~]=kde(R3(:, 1), den(4,1), 0, 4); plot(xmesh3, density3, 'Color',colour_green, 'DisplayName', names{3});
xlabel('\theta_1');


subplot(1, 3, 2); hold on; grid on;
[~,density,xmesh,~]=kde(S_rej(:, 2), den(1,2), 0, 0.04); plot(xmesh, density, 'LineWidth', 1, 'Color', rgb('Black'), 'DisplayName', 'Reference');
[~,density1,xmesh1,~]=kde(R1(:, 2), den(2,2), 0, 0.04); plot(xmesh1, density1, 'Color',colour_darkblue, 'DisplayName', names{1});
[~,density2,xmesh2,~]=kde(R2(:, 2), den(3,2), 0, 0.04); plot(xmesh2, density2, 'Color',colour_peach, 'DisplayName', names{2});
[~,density3,xmesh3,~]=kde(R3(:, 2), den(4,2), 0, 0.04); plot(xmesh3, density3, 'Color',colour_green, 'DisplayName', names{3});
xlabel('\theta_2');

subplot(1, 3, 3);  hold on; grid on;
[~,density,xmesh,~]=kde(S_rej(:, 3), den(1,3), 0, 4); plot(xmesh, density, 'LineWidth', 1, 'Color', rgb('Black'), 'DisplayName', 'Reference');
[~,density1,xmesh1,~]=kde(R1(:, 3), den(2,3), 0, 4); plot(xmesh1, density1, 'Color',colour_darkblue, 'DisplayName', names{1});
[~,density2,xmesh2,~]=kde(R2(:, 3), den(3,3), 0, 4); plot(xmesh2, density2, 'Color',colour_peach, 'DisplayName', names{2});
[~,density3,xmesh3,~]=kde(R3(:, 3), den(4,3), 0, 4); plot(xmesh3, density3, 'Color',colour_green, 'DisplayName', names{3});
xlabel('\theta_3');


% xlabel('lag'); ylabel('sample autocorrelation')
%title(titl, 'FontWeight', 'Normal');
lgd = legend('show');
lgd.FontSize=7;
%lgd.Orientation='horizontal';
lgd.Location='north east';
end
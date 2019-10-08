function plot_three_means(R1, R2, R3, S_rej, names)
% Plot sample acf of three algorithm outputs for comparison
% Algorithms:   1) Multiple processor APTMC
%               2) Single processor APTMC
%               3) Standard AMC

if ~exist('names', 'var')
    names = {'ABC', 'ABC-PTMC-1', 'ABC-APTMC-1'};
end

colour_darkblue = [1 17 181] ./ 255; %ABC-PTMC-1
colour_peach = [251 111 66] ./ 255; %ABC-APTMC-1
colour_green = [12 195 82] ./ 255;  %ABC

% Get posterior means by iteration
n1=size(R1,1); 
mns1 = cumsum(R1)./[1:n1; 1:n1; 1:n1]';

n2=size(R2,1); 
mns2 = cumsum(R2)./[1:n2; 1:n2; 1:n2]';

n3=size(R3,1);
mns3 = cumsum(R3)./[1:n3; 1:n3; 1:n3]';

rej_mean = mean(S_rej); 

%THETA_1
figure; 
subplot(3, 1, 1); hold on;
plot(mns1(:, 1), 'Color',colour_green, 'DisplayName', names{1}); 
plot(mns2(:, 1), 'Color',colour_darkblue, 'DisplayName', names{2})
plot(mns3(:, 1), 'Color',colour_peach, 'DisplayName', names{3})
ylabel('mean \theta_1'); xlim([0 n2]);
hline=refline(0, rej_mean(1)); hline.DisplayName = 'reference'; hline.Color =rgb('Black');
xlabel('samples');
lgd = legend('show');
lgd.FontSize=7;
lgd.Location='north east';

% THETA_2
subplot(3, 1, 2); hold on;
plot(mns1(:, 2), 'Color',colour_green, 'DisplayName', names{1}); 
plot(mns2(:, 2), 'Color',colour_darkblue, 'DisplayName', names{2})
plot(mns3(:, 2), 'Color',colour_peach, 'DisplayName', names{3})
xlim([0 n2]);
ylabel('mean \theta_2');
hline=refline(0, rej_mean(2)); hline.Color = rgb('Black');
xlabel('samples');

%THETA_3
subplot(3, 1, 3); hold on
plot(mns1(:, 3), 'Color',colour_green, 'DisplayName', names{1}); 
plot(mns2(:, 3), 'Color',colour_darkblue, 'DisplayName', names{2})
plot(mns3(:, 3), 'Color',colour_peach, 'DisplayName', names{3});
ylabel('mean \theta_3'); xlim([0 n2]);
hline=refline(0, rej_mean(3)); hline.Color = rgb('Black');

xlabel('samples');


end
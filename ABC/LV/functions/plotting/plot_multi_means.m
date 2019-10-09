function plot_multi_means(R, S_rej)
% Plot sample acf of three algorithm outputs for comparison
% Algorithms:   1) Multiple processor APTMC
%               2) Single processor APTMC
%               3) Standard AMC


ntimes = max(size(R)); mns = cell(1, ntimes);
colour_darkblue = [1 17 181] ./ 255; %ABC-PTMC-1
colour_peach = [251 111 66] ./ 255; %ABC-APTMC-1
colour_green = [12 195 82] ./ 255;  %ABC

% Get posterior means by iteration
for idx=1:ntimes
    n1=size(R{idx},1);
    mns{idx} = cumsum(R{idx})./[1:n1; 1:n1; 1:n1]';
end

rej_mean = mean(S_rej); 

%THETA_1
figure; 
subplot(3, 1, 1); hold on;
for idx=1:ntimes
    plot(mns{idx}(:, 1), 'Color', colour_green);
end
ylabel('mean \theta_1'); %xlim([0 n2]);
hline=refline(0, rej_mean(1)); hline.DisplayName = 'reference'; hline.Color =rgb('Black');
xlabel('samples');
% lgd = legend('show');
% lgd.FontSize=7;
% lgd.Location='north east';

% THETA_2
subplot(3, 1, 2); hold on;
for idx=1:ntimes
    plot(mns{idx}(:, 2), 'Color', colour_darkblue);
end
ylabel('mean \theta_2'); % xlim([0 n2]);
hline=refline(0, rej_mean(2)); hline.Color = rgb('Black');
xlabel('samples');

%THETA_3
subplot(3, 1, 3); hold on
for idx=1:ntimes
    plot(mns{idx}(:, 3), 'Color', colour_peach);
end

ylabel('mean \theta_3'); %xlim([0 n2]);
hline=refline(0, rej_mean(3)); hline.Color = rgb('Black');

xlabel('samples');


end
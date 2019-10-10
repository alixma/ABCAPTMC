function plot_essovertime(W, R, iat, T, cis)


TMA = []; TMS = []; TMM=[];
for w=1:W
    % ABC-APTMC-1
    TMA = [TMA; R{1}{w}(:, 3:4)];
    % ABC-PTMC-1
    TMS = [TMS; R{2}{w}(:, 3:4)];
    % Standard ABC
    TMM = [TMM; R{3}{w}(R{3}{w}(:, 1)>0, :)];
end

col = [[12 195 82] ./ 255; rgb('DarkBlue');  rgb('Red')];
scol = [rgb('DarkSeaGreen'); rgb('LightBlue');  rgb('Goldenrod')];
sz = .1; alph = .5; lw = 1.5;
figure;
for i=1:3
    p_a = fitlm(TMA(:, 1), TMA(:, 2)/iat.anytime(i), 'Intercept', false); [fit_a, int_a] = predict(p_a, TMA(:, 1));
    p_e = fitlm(TMS(:, 1), TMS(:, 2)/iat.stdexchange(i), 'Intercept', false); [fit_e, int_e] = predict(p_e, TMS(:, 1));
    p_m = fitlm(TMM(:, 1), TMM(:, 2)/iat.standard(i), 'Intercept', false); [fit_m, int_m] = predict(p_m, TMM(:, 1));
    subplot(3, 1, i);
    hold on;
    % ABC-APTMC-1
    g1 = scatter(TMA(:, 1), TMA(:, 2)/iat.anytime(i), sz, 'Marker', '.',...
        'MarkerEdgeColor', scol(3,:),...
        'MarkerFaceColor', scol(3,:),...
        'MarkerFaceAlpha', alph/5,...
        'MarkerEdgeAlpha', alph/5,...
        'DisplayName','ABC-APTMC');
    % ABC-PTMC-1
    g2 = scatter(TMS(:, 1), TMS(:, 2)/iat.stdexchange(i), sz, 'MarkerEdgeColor', scol(2,:),...
        'MarkerFaceColor', scol(2,:),...
        'MarkerFaceAlpha', alph,...
        'MarkerEdgeAlpha', alph,...
        'DisplayName','ABC-PTMC');
    % Standard ABC
    g3 = scatter(TMM(:, 1), TMM(:, 2)/iat.standard(i), sz, 'Marker', 'd',...
        'MarkerEdgeColor', scol(1,:),...
        'MarkerFaceColor', scol(1,:),...
        'MarkerFaceAlpha', alph,...
        'MarkerEdgeAlpha', alph,...
        'DisplayName','Standard');
    h1 = plot(TMA(:, 1), fit_a, 'Color', col(3,:), 'LineWidth', lw, 'DisplayName','ABC-APTMC fit');
    if cis
        h1_1 = plot(TMA(:, 1), int_a(:, 1), 'LineStyle', '--', 'Color', col(3,:), 'DisplayName','ABC-APTMC CI');
        h1_2 = plot(TMA(:, 1), int_a(:, 2), 'LineStyle', '--', 'Color', col(3,:));
    end
    h2 = plot(TMS(:, 1), fit_e, 'Color', col(2,:), 'LineWidth', lw, 'DisplayName','ABC-PTMC fit');
    if cis
        h2_1 = plot(TMS(:, 1), int_e(:, 1), 'LineStyle', '--', 'Color', col(2,:), 'DisplayName','ABC-PTMC CI');
        h2_2 = plot(TMS(:, 1), int_e(:, 2), 'LineStyle', '--', 'Color', col(2,:));
    end
    h3 = plot(TMM(:, 1), fit_m, 'Color', col(1,:), 'LineWidth', lw, 'DisplayName','Standard fit');
    if cis
        h3_1 = plot(TMM(:, 1), int_m(:, 1), 'LineStyle', '--', 'Color', col(1,:), 'DisplayName','Standard CI');
        h3_2 = plot(TMM(:, 1), int_m(:, 2), 'LineStyle', '--', 'Color', col(1,:));
    end
    xlabel('time'); xlim([0 T]);
    ylabel('ESS')
    title(sprintf('ESS over time, \\theta_%d', i))
    if cis
        legend([h1 h1_1 h2 h2_1 h3 h3_1], 'Location', 'northwest');
    else
        legend('Location', 'northwest');
    end
end
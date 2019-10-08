function plot_allchains(S_in, params, observations, bins, nden, ax)
% Plot all chains output from the ABC-APTMC algorithm

if(~exist('ax', 'var'))
    ax=[-7.5 7.5 0.0 0.5];
end
if(~exist('bins', 'var'))
    bins=30;
end
figure
np=2^10;
xx=linspace(ax(1), ax(2), np);
xxden = linspace(ax(1), ax(2), nden);
yy = normpdf(xx, params.posterior(1), params.posterior(2));
C = size(S_in, 2);


for k=1:params.K
    h = subplot(2, round(params.K/2), k);
    pos = get(h, 'pos');
    pos(1) = pos(1) - 0.04;
    pos(2) = pos(2) - 0.05;
    pos(3) = pos(3) + 0.04;
    pos(4) = pos(4) + 0.04;
    
    set(h, 'pos', pos);    
    if C==params.K
        S0 = S_in{k};
    else
        S0=[];
        for c=1:C
            S0 = [S0; S_in{c}{k}];
        end
    end
    %numerical integration to obtain tempered distributions
    F = @(x) likelihood(observations.y, x, params.epsilon(k), observations.sigma).*normpdf(x, params.prior(1), params.prior(2));
    Z = integral(F, ax(1), ax(2));
    zz = 1/Z*F(xx);
    [~, den, ~, ~] = kde(S0, nden, ax(1), ax(2));
    
    plot(xx, yy,...
    'Color', 'green',...
    'LineWidth',1.4);
    hold on

    h=histogram(S0, 'Normalization','pdf');
    h.NumBins = bins;
    set(h,'FaceColor',rgb('LightBlue'),'EdgeColor',rgb('DodgerBlue'));
    
    plot(xx, zz,...
    'Color', 'red',...
    'LineWidth',1.2);
    plot(xxden, den,...
    'Color', rgb('BlueViolet'),...
     'linestyle', '--',...
    'LineWidth',1.2);
    grid on
    axis(ax);
    title([sprintf('Posterior for chain #%d, radius ', k), char(949), sprintf(' = %f', params.epsilon(k))], 'FontWeight', 'Normal');
    if(k~=1 && k~=round(params.K/2)+1)
        set(gca, 'YTickLabel', {});
    end
    if(k<params.K/2+1)
        set(gca, 'XTickLabel', {});
    end
end
end

function plot_allchains_both(S_in, params, observations, nden, col, ax, names)
% Plot all chains output from the ABC-APTMC algorithm
% S_in must already have combined all parallel runs

if(~exist('ax', 'var'))
    ax=[-7.5 7.5 0.0 0.5];
end
if(~exist('col', 'var'))
    col = [rgb('Lime'); rgb('Tomato')];
end
if(~exist('names', 'var'))
    names={'corrected', 'uncorrected'};
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
    %numerical integration to obtain tempered distributions
    F = @(x) likelihood(observations.y, x, params.epsilon(k), observations.sigma).*normpdf(x, params.prior(1), params.prior(2));
    Z = integral(F, ax(1), ax(2));
    zz = 1/Z*F(xx);
    plot(xx, yy,...
    'Color', rgb('LightGray'),...
    'LineWidth',1.5,...
    'DisplayName','cold posterior');
    hold on
     plot(xx, zz,...
    'Color', rgb('DimGray'),...
    'LineWidth', 2,...
    'DisplayName','chain posterior');
    for c=1:C
        S0 = S_in{c}{k};
        [~, den, ~, ~] = kde(S0, nden, ax(1), ax(2));
         plot(xxden, den,...
            'Color', col(c,:),...
             'linestyle', '-.',...
            'LineWidth', 1.2,...
            'DisplayName',names{c});
    end        
  
    grid on
    axis(ax);
    
    title([sprintf('Posterior for chain #%d, radius ', k), char(949), sprintf(' = %f', params.epsilon(k))], 'FontWeight', 'Normal');
    if(k~=1 && k~=round(params.K/2)+1)
        set(gca, 'YTickLabel', {});        
    end
    if k==round(params.K/4)
        lgd = legend;
        newPosition = [0.448281250737296,0.475020555474654,0.102604165192072,0.037512845885227];
        newUnits = 'normalized';
        set(lgd,'Position', newPosition,'Units', newUnits, 'Orientation', 'vertical', 'NumColumns', 2);
    end
    if(k<params.K/2+1)
        set(gca, 'XTickLabel', {});
    end
end
end

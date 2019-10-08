function MA_plot_allchains(Theta, n, ne, theta_true, params, Y, colour, sz, ax)

if(~exist('ax', 'var'))
    ax=[-2 2 -1 1];
end

if(~exist('sz', 'var'))
    sz=10;
end

figure
xx = [-2 0 2];
yy = [1 -1 1];

for k=1:params.K
    fprintf('Plotting chain %d \n', k)
    h = subplot(2, round(params.K/2), k);
    
    h = fill(xx,yy,rgb('Lavender')); hold on;
    
    S1 = Theta((ne(k)+1):n(k), 1, k);
    S2 = Theta((ne(k)+1):n(k), 2, k);
    p=gkde2([S1,S2]);
    
    [cx, cy] = fixxy(p.x, p.y);
      
    S = MA_likelihood2(cx, cy, Y, 1.749528420996907e-25);
    
    
    scatter(S1, S2, sz, colour, 'filled');
    %contour(cx, cy, S, 'LineColor', 'black', 'LineStyle', '-.');    
    zlevs = get_levels(S, 5);
    p.zlevs = get_levels(p.pdf, 5);
    contour(cx, cy, S, 'LineColor', 'black', 'LevelList', zlevs, 'LineStyle', '-.', 'LineWidth', .1);
    contour(p.x,p.y,p.pdf, 'LevelList', p.zlevs);
    scatter(theta_true(1), theta_true(2), 50, rgb('LawnGreen'), 'filled', 'pentagram');
    
    axis(ax);
    title([sprintf('Chain #%d, ', k), char(949), sprintf(' = %f', params.epsilon(k))], 'FontWeight', 'Normal');
    %     if(k~=1 && k~=round(params.K/2)+1)
    %         set(gca, 'YTickLabel', {});
    %     else
    ylabel('\theta_2');
    %     end
    %     if(k<params.K/2+1)
    %         set(gca, 'XTickLabel', {});
    %     else
    xlabel('\theta_1');
    %     end
end

end
function  subfun_plot_scatter_regression( x , y , xlim , ylim )

% figure; hold on;
% draw


mdl = fitlm(x, y);
xfit = linspace(min(x), max(x), 200)';
[yhat, yCI] = predict(mdl, xfit);  

fill([xfit; flipud(xfit)], [yCI(:,1); flipud(yCI(:,2))], ...
     [0.8 0.8 0.8], 'EdgeColor', 'none');       % 灰色置信区间
scatter(x, y, 10, [0.7 0.7 0.7], 'filled');
plot(xfit, yhat, '-', 'color', [0 0 0], 'LineWidth', 1);

set(gca,'YDir','normal'); 
set(gca, 'ylim', ylim, 'yTick', ylim,...
         'xlim', xlim, 'xTick', xlim,...
         'TickDir', 'out');
axis square


end


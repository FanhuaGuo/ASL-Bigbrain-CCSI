function  subfun_plot_scatter_density_isopotential( x , y , valueorrank , ifdensity , ifisopotential , ifregressline )
%   its a subfunction to plot the correlation results with scatter, density
% and the isopotential line.

if valueorrank == 2
    [~,xr] = sort(x);
    [~,yr] = sort(y);
    y = zeros(numel(yr),1);
    x = zeros(numel(xr),1);
    for i = 1:numel(xr)
        x(xr(i)) = i;
        y(yr(i)) = i;
    end
end

% grid point
[xq, yq] = meshgrid(linspace(min(x),max(x),200), linspace(min(y),max(y),200));

% kernel density estimation
f = ksdensity([x,y],[xq(:),yq(:)]);
f = reshape(f, size(xq));

threshold = [0 max(f(:))];
Gcolors = [linspace(1,0.5,1000)' linspace(1,0.5,1000)' linspace(1,0.5,1000)'];

% figure; hold on;
% draw
if ifdensity
    imagesc(linspace(min(x),max(x),200), linspace(min(y),max(y),200), f, threshold);
end
if ifisopotential
    contour(xq, yq, f, 15, 'k'); 
end
scatter(x, y, 10, [0 0 0], 'filled');
if ifregressline
    [b] = regress(y, [ones(size(x)) x]); % 回归
    x_line = [min(x), max(x)];
    y_line = b(1) + b(2)*x_line;
    plot(x_line, y_line, '-', 'color', [0 0 0], 'LineWidth', 3);
end
colormap(Gcolors); 
set(gca,'YDir','normal'); 
set(gca, 'ylim', [min(y),max(y)], 'yTick', [],...
         'xlim', [min(x),max(x)], 'xTick', []);
set(gca, 'Box', 'off');
ax = gca;
ax.XColor = 'none'; 
ax.YColor = 'none'; 
ax.TickDir = 'out';
axis square


end


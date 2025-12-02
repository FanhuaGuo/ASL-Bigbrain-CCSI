function  subfun_plotting_corr_Pcorr_heatmaps( draw_r , draw_p , draw_label , Threshold , ifplotcommon )
%% draw corr with maps (like gradient)
% disposal draw-data
xticklabel = {'common','L6','L5','L4','L3','L2','L1'};
if ~ifplotcommon
    draw_r(:,3) = [];
    draw_p(:,3) = [];
    xticklabel = {'L6','L5','L4','L3','L2','L1'};
end

% colormap
n = 256;
r = [linspace(0.3,1,n/2) linspace(1,0.65,n/2)];
g = [linspace(0.3,1,n/2) linspace(1,0.1,n/2)];
b = [linspace(0.6,1,n/2) linspace(1,0.15,n/2)];
cmap = [r(:) g(:) b(:)];
cmap(128,:) = [];

% draw CBF-score and CCSI
figure('Color',[1 1 1],'Position',[0 0 300 numel(draw_label)*50],'Units','pixels');
hold on;
Matrix_draw = draw_r(:,1:2);
imagesc(Matrix_draw,Threshold);
colormap(cmap);
colorbar
set(gca,'xLim',[0.5 2+0.5],'xTick',[1 2],'xTickLabel',{'CBF-score','CCSI'},'FontName', 'Arial');
set(gca,'yLim',[0.5 numel(draw_label)+0.5],'yTick',[1:numel(draw_label)],'yTickLabel',draw_label,'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.XColor = 'k'; % 保持刻度数字
ax.YColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';
% 删除 axis 的 Label
xlabel('cortical depth');
ylabel('normalized value');
% 隐藏多余的轴元素
ax.XRuler.Axle.Visible = 'off';
ax.YRuler.Axle.Visible = 'off';

for i = 1:size(Matrix_draw,1)
    for j = 1:2
        if draw_p(i,j) < 0.001
            txt = '<0.001';
        else
            txt = sprintf('%.3f', draw_p(i,j));
        end
        text(j, i, txt, 'HorizontalAlignment', 'center', ...
            'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
    end
end
% 添加显著性边框
for i = 1:size(Matrix_draw,1)
    for j = 1:2
        if draw_p(i,j) < 0.05
            rectangle('Position', [j-0.5 i-0.5 1 1], 'EdgeColor', 'k', 'LineWidth', 1.5);
        end
    end
end 


% draw laminar CBF-score
figure('Color',[1 1 1],'Position',[0 0 100*numel(xticklabel) numel(draw_label)*50],'Units','pixels');
hold on;
Matrix_draw = draw_r(:,3:end);
imagesc(Matrix_draw,Threshold);
colormap(cmap);
colorbar
set(gca,'xLim',[0.5 numel(xticklabel)+0.5],'xTick',[1:numel(xticklabel)],'xTickLabel',xticklabel, 'FontName', 'Arial');
set(gca,'yLim',[0.5 numel(draw_label)+0.5],'yTick',[]);
set(gca, 'Box', 'off');
ax = gca;
ax.XColor = 'k'; % 保持刻度数字
ax.YColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';
% 删除 axis 的 Label
xlabel('cortical depth');
ylabel('normalized value');
% 隐藏多余的轴元素
ax.XRuler.Axle.Visible = 'off';
ax.YRuler.Axle.Visible = 'off';
for i = 1:size(Matrix_draw,1)
    for j = 1:numel(xticklabel)
        if draw_p(i,j+2) < 0.001
            txt = '<0.001';
        else
            txt = sprintf('%.3f', draw_p(i,j+2));
        end
        text(j, i, txt, 'HorizontalAlignment', 'center', ...
            'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
    end
end
for i = 1:size(Matrix_draw,1)
    for j = 1:numel(xticklabel)
        if draw_p(i,j+2) < 0.05
            rectangle('Position', [j-0.5 i-0.5 1 1], 'EdgeColor', 'k', 'LineWidth', 1.5);
        end
    end
end 

end
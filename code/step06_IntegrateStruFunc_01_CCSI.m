clc
clear all
close all
set(0,'defaultfigurecolor',[1 1 1]);% set figure background color
if ~contains(path, '/Users/guofanhua/Desktop/gfh/tools/MatlabTools/BrainSpace-latest/matlab')
    addpath(genpath('/Users/guofanhua/Desktop/gfh/tools/MatlabTools/BrainSpace-latest/matlab'));
end

%% ======================= compare regional CCSI

%% input data
step01_CBFandCCSImaps_00_init;
step01_CBFandCCSImaps_00_ComputeBasicResults

ccsi = load([Group_dir 'CBF_Mean_PScs_L12_bil.1D']);
ccsi_null = load([Group_dir 'CBF_Mean_PScs_L12_bil_variogram.1D']);
cbf = load([Group_dir 'CBF_PC1_Val_L12_bil.1D']);
cbf_null = load([Group_dir 'CBF_PC1_Val_L12_bil_variogram.1D']);

bb_g = BB_gradient;
fc1_g = all_other_maps.maps{38};
allp = zeros(9,1);


% write to SourceData Fig8abd
outSD = [bb_g fc1_g bb_g];
outSD(outSD(:,3)<0,3) = 0;
outSD(outSD(:,3)>0,3) = 1;


%% whole bb gradient
% close all
% clc
mask = CBFsnr_mean_bil>4;  % bb_g>0 

% plotting
figure('Color',[1 1 1],'Position',[0 0 300 300],'Units','pixels');
hold on;
% subfun_plot_scatter_density_isopotential(stru,func, 2 , 1 , 1 , 0 )
subfun_plot_scatter_regression(bb_g(mask),fc1_g(mask), [-0.12 0.12], [-8 8]);


% set
func = zscore(fc1_g(mask)); 
stru = zscore(bb_g(mask)); 
meta_cbf = zscore(cbf(mask)); 
null1 = zscore(cbf_null(mask,:));
meta_ccsi = zscore(ccsi(mask)); 
null2 = zscore(ccsi_null(mask,:));
fprintf('For whole MPC gradient:\n\n');


% VIF
% 计算相关系数矩阵的逆，其对角线元素即为 VIF
% 或者使用统计工具箱自带函数
R = corrcoef([stru, meta_ccsi]);
vif_values = diag(inv(R))'; 
fprintf('VIF values: G_MPC = %.2f, CCSI = %.2f\n', vif_values(1), vif_values(2));




% ---------- Step 1: baseline model ----------
X1 = [ones(length(stru),1), stru];   % 截距 + Stru
[b1,~,~,~,stats1] = regress(func, X1);
R2_1 = stats1(1); % R^2
fprintf('Basic model:\nStep1: Func = %.3f + %.3f*Stru, R^2=%.3f, p=%.3g\n', b1(1), b1(2), R2_1, stats1(3));
residuals_base = func - b1(1) - b1(2)*stru;


% % 第二步：可视化残差与 CCSI 的关系
% figure;
% scatter(meta_ccsi, residuals_base, 'filled', 'MarkerFaceAlpha', 0.5);
% hold on;
% lsline; % 添加拟合线
% xlabel('CCSI (Perfusion-Structure Coupling)');
% ylabel('Residuals from G_{MPC} model');
% title('Does CCSI explain the "unexplained" functional variance?');
% 
% % 第三步：量化 CCSI 对残差的解释力
% mdl_resid = fitlm(meta_ccsi, residuals_base);
% [r_val, p_val] = corr(meta_ccsi, residuals_base, 'type', 'Pearson');
% 
% fprintf('CCSI vs Residuals: r = %.3f, p = %.4f\n', r_val, p_val);




% ---------- Step 2: additional model (ccsi) ----------
X2 = [ones(length(stru),1), stru, meta_ccsi]; % 截距 + Stru + Meta
[b2,~,~,~,stats2] = regress(func, X2);
R2_2 = stats2(1); % R^2
fprintf('\nAdd CCSI:\nStep3: Func = %.3f + %.3f*Stru + %.3f*CCSI, R^2=%.3f, p=%.3g\n', ...
    b2(1), b2(2), b2(3), R2_2, stats2(3));

% spin-null test
nperm = size(null2,2);
deltaR2_null = zeros(nperm,1);
deltab_null = zeros(nperm,1);
bm_null = zeros(nperm,1);
for ip = 1:nperm
    X2_perm = [ones(length(stru),1), stru, null2(:,ip)];
    [bn,~,~,~,stats_perm] = regress(func, X2_perm);
    deltaR2_null(ip) = stats_perm(1) - R2_1;
    deltab_null(ip) = bn(2) - b1(2);
    bm_null(ip) = bn(3);
end

deltaR2 = R2_2 - R2_1; 
p_permr2 = mean(deltaR2_null >= deltaR2);
fprintf('ΔR² (Increased explanatory power) = %.3f,  ', deltaR2);
fprintf('ΔR² spin-null p-value = %.3g\n', p_permr2);

deltab = b2(2) - b1(2);
p_permb = mean(abs(deltab_null) >= abs(deltab));
fprintf('Δbeta (increased stru-func coupling) = %.3f,  ', deltab);
fprintf('Δbeta spin-null p-value = %.3g\n', p_permb);

allp(1) = p_permr2;
allp(2) = p_permb;

bm = b2(3);
p_permb = mean(abs(bm_null) >= abs(bm));
fprintf('ccsi beta = %.3f,  ', bm);
fprintf('ccsi beta spin-null p-value = %.3g\n', p_permb);

allp(7) = p_permb;



% plot
figure('Color',[1 1 1],'Position',[0 0 600 300],'Units','pixels');
subplot(141); hold on;
bar(1, R2_1, 'FaceColor', 'none', 'EdgeColor', [0, 200, 180]/255,'linewidth',1);
bar(3, R2_2, 'FaceColor', 'none', 'EdgeColor', [150, 103, 185]/255,'linewidth',1);
set(gca,'xLim',[0 4],'xTick',[1 3], 'xTickLabel', {'basic R2','add CCSI R2'},'FontName', 'Arial');
set(gca,'yLim',[0.2 0.3],'yTick',[0.2 0.3], 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';

subplot(142); hold on;
% swarmchart([ones(numel(deltaR2_null),1)],deltaR2_null,5,[0.5 0.5 0.5],'filled','XJitter','rand','XJitterWidth',0.5);
swarmchart([ones(numel(deltaR2_null),1)],deltaR2_null,5,[0.5 0.5 0.5],'filled','XJitterWidth',0.5);
scatter(1,deltaR2,20,[1 0 0],'filled')
set(gca,'xLim',[0.5 1.5],'xTick',[1], 'xTickLabel', {'ΔR2 null'},'FontName', 'Arial');
set(gca,'yLim',[0 0.25],'yTick',[0 0.25], 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';

subplot(143); hold on;
bar(1, b1(2), 'FaceColor', 'none', 'EdgeColor', [0, 200, 180]/255,'linewidth',1);
bar(3, b2(2), 'FaceColor', 'none', 'EdgeColor', [150, 103, 185]/255,'linewidth',1);
set(gca,'xLim',[0 4],'xTick',[1 3], 'xTickLabel', {'basic beta','add CCSI beta'},'FontName', 'Arial');
set(gca,'yLim',[0.4 0.6],'yTick',[0.4 0.6], 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';

subplot(144); hold on;
% swarmchart([ones(numel(deltab_null),1)],abs(deltab_null),5,[0.5 0.5 0.5],'filled','XJitter','rand','XJitterWidth',0.5);
swarmchart([ones(numel(deltab_null),1)],abs(deltab_null),5,[0.5 0.5 0.5],'filled','XJitterWidth',0.5);
scatter(1,deltab,20,[1 0 0],'filled')
set(gca,'xLim',[0.5 1.5],'xTick',[1], 'xTickLabel', {'Δbeta null'},'FontName', 'Arial');
set(gca,'yLim',[0 0.2],'yTick',[0 0.2], 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';






% plot 20260406
figure('Color',[1 1 1],'Position',[0 0 600 300],'Units','pixels');
subplot(131); hold on;
% swarmchart([ones(numel(deltaR2_null),1)],deltaR2_null,5,[0.5 0.5 0.5],'filled','XJitter','rand','XJitterWidth',0.5);
swarmchart([ones(numel(deltaR2_null),1)],deltaR2_null,5,[0.5 0.5 0.5],'filled','XJitterWidth',0.5);
scatter(1,deltaR2,20,[1 0 0],'filled')
set(gca,'xLim',[0.5 1.5],'xTick',[1], 'xTickLabel', {'ΔR2 null'},'FontName', 'Arial');
set(gca,'yLim',[0 0.25],'yTick',[0 0.25], 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';

subplot(132); hold on;
% swarmchart([ones(numel(deltaR2_null),1)],deltaR2_null,5,[0.5 0.5 0.5],'filled','XJitter','rand','XJitterWidth',0.5);
swarmchart([ones(numel(deltab_null),1)],abs(deltab_null),5,[0.5 0.5 0.5],'filled','XJitterWidth',0.5);
scatter(1,abs(deltab),20,[1 0 0],'filled')
set(gca,'xLim',[0.5 1.5],'xTick',[1], 'xTickLabel', {'Δbeta null'},'FontName', 'Arial');
set(gca,'yLim',[0 0.2],'yTick',[0 0.2], 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';

subplot(133); hold on;
% swarmchart([ones(numel(deltaR2_null),1)],deltaR2_null,5,[0.5 0.5 0.5],'filled','XJitter','rand','XJitterWidth',0.5);
swarmchart([ones(numel(bm_null),1)],abs(bm_null),5,[0.5 0.5 0.5],'filled','XJitterWidth',0.5);
scatter(1,abs(bm),20,[1 0 0],'filled')
set(gca,'xLim',[0.5 1.5],'xTick',[1], 'xTickLabel', {'ccsi beta null'},'FontName', 'Arial');
set(gca,'yLim',[0 0.6],'yTick',[0 0.6], 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';



%% Primary cortex: bb<0
% close all
% clc
mask = CBFsnr_mean_bil>4  &  bb_g<0 ;  

% plotting
figure('Color',[1 1 1],'Position',[0 0 300 300],'Units','pixels');
hold on;
% subfun_plot_scatter_density_isopotential(stru,func, 2 , 1 , 1 , 0 )
subfun_plot_scatter_regression(bb_g(mask),fc1_g(mask), [-0.12 0.12], [-8 8]);


% set
func = zscore(fc1_g(mask)); 
stru = zscore(bb_g(mask)); 
meta_cbf = zscore(cbf(mask)); 
null1 = zscore(cbf_null(mask,:));
meta_ccsi = zscore(ccsi(mask)); 
null2 = zscore(ccsi_null(mask,:));
fprintf('For bb gradient<0:\n\n');


% VIF
% 计算相关系数矩阵的逆，其对角线元素即为 VIF
% 或者使用统计工具箱自带函数
R = corrcoef([stru, meta_ccsi]);
vif_values = diag(inv(R))'; 
fprintf('VIF values: G_MPC = %.2f, CCSI = %.2f\n', vif_values(1), vif_values(2));




% ---------- Step 1: baseline model ----------
X1 = [ones(length(stru),1), stru];   % 截距 + Stru
[b1,~,~,~,stats1] = regress(func, X1);
R2_1 = stats1(1); % R^2
fprintf('Basic model:\nStep1: Func = %.3f + %.3f*Stru, R^2=%.3f, p=%.3g\n', b1(1), b1(2), R2_1, stats1(3));





% ---------- Step 2: additional model (ccsi) ----------
X2 = [ones(length(stru),1), stru, meta_ccsi]; % 截距 + Stru + Meta
[b2,~,~,~,stats2] = regress(func, X2);
R2_2 = stats2(1); % R^2
fprintf('\nAdd CCSI:\nStep3: Func = %.3f + %.3f*Stru + %.3f*CCSI, R^2=%.3f, p=%.3g\n', ...
    b2(1), b2(2), b2(3), R2_2, stats2(3));

% spin-null test
nperm = size(null2,2);
deltaR2_null = zeros(nperm,1);
deltab_null = zeros(nperm,1);
bm_null = zeros(nperm,1);
for ip = 1:nperm
    X2_perm = [ones(length(stru),1), stru, null2(:,ip)];
    [bn,~,~,~,stats_perm] = regress(func, X2_perm);
    deltaR2_null(ip) = stats_perm(1) - R2_1;
    deltab_null(ip) = bn(2) - b1(2);
    bm_null(ip) = bn(3);
end

deltaR2 = R2_2 - R2_1; 
p_permr2 = mean(deltaR2_null >= deltaR2);
fprintf('ΔR² (Increased explanatory power) = %.3f,  ', deltaR2);
fprintf('ΔR² spin-null p-value = %.3g\n', p_permr2);

deltab = b2(2) - b1(2);
p_permb = mean(abs(deltab_null) >= abs(deltab));
fprintf('Δbeta (increased stru-func coupling) = %.3f,  ', deltab);
fprintf('Δbeta spin-null p-value = %.3g\n', p_permb);

allp(3) = p_permr2;
allp(4) = p_permb;

bm = b2(3);
p_permb = mean(abs(bm_null) >= abs(bm));
fprintf('ccsi beta = %.3f,  ', bm);
fprintf('ccsi beta spin-null p-value = %.3g\n', p_permb);

allp(8) = p_permb;



% plot
figure('Color',[1 1 1],'Position',[0 0 600 300],'Units','pixels');
subplot(141); hold on;
bar(1, R2_1, 'FaceColor', 'none', 'EdgeColor', [0, 200, 180]/255,'linewidth',1);
bar(3, R2_2, 'FaceColor', 'none', 'EdgeColor', [150, 103, 185]/255,'linewidth',1);
set(gca,'xLim',[0 4],'xTick',[1 3], 'xTickLabel', {'basic R2','add CCSI R2'},'FontName', 'Arial');
set(gca,'yLim',[0 0.2],'yTick',[0 0.2], 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';

subplot(142); hold on;
% swarmchart([ones(numel(deltaR2_null),1)],deltaR2_null,5,[0.5 0.5 0.5],'filled','XJitter','rand','XJitterWidth',0.5);
swarmchart([ones(numel(deltaR2_null),1)],deltaR2_null,5,[0.5 0.5 0.5],'filled','XJitterWidth',0.5);
scatter(1,deltaR2,20,[1 0 0],'filled')
set(gca,'xLim',[0.5 1.5],'xTick',[1], 'xTickLabel', {'ΔR2 null'},'FontName', 'Arial');
set(gca,'yLim',[0 0.4],'yTick',[0 0.4], 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';

subplot(143); hold on;
bar(1, b1(2), 'FaceColor', 'none', 'EdgeColor', [0, 200, 180]/255,'linewidth',1);
bar(3, b2(2), 'FaceColor', 'none', 'EdgeColor', [150, 103, 185]/255,'linewidth',1);
set(gca,'xLim',[0 4],'xTick',[1 3], 'xTickLabel', {'basic beta','add CCSI beta'},'FontName', 'Arial');
set(gca,'yLim',[0.3 0.4],'yTick',[0.3 0.4], 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';

subplot(144); hold on;
% swarmchart([ones(numel(deltab_null),1)],abs(deltab_null),5,[0.5 0.5 0.5],'filled','XJitter','rand','XJitterWidth',0.5);
swarmchart([ones(numel(deltab_null),1)],abs(deltab_null),5,[0.5 0.5 0.5],'filled','XJitterWidth',0.5);
scatter(1,deltab,20,[1 0 0],'filled')
set(gca,'xLim',[0.5 1.5],'xTick',[1], 'xTickLabel', {'Δbeta null'},'FontName', 'Arial');
set(gca,'yLim',[0 0.3],'yTick',[0 0.3], 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';








% plot 20260406
figure('Color',[1 1 1],'Position',[0 0 600 300],'Units','pixels');
subplot(131); hold on;
% swarmchart([ones(numel(deltaR2_null),1)],deltaR2_null,5,[0.5 0.5 0.5],'filled','XJitter','rand','XJitterWidth',0.5);
swarmchart([ones(numel(deltaR2_null),1)],deltaR2_null,5,[0.5 0.5 0.5],'filled','XJitterWidth',0.5);
scatter(1,deltaR2,20,[1 0 0],'filled')
set(gca,'xLim',[0.5 1.5],'xTick',[1], 'xTickLabel', {'ΔR2 null'},'FontName', 'Arial');
set(gca,'yLim',[0 0.4],'yTick',[0 0.4], 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';

subplot(132); hold on;
% swarmchart([ones(numel(deltaR2_null),1)],deltaR2_null,5,[0.5 0.5 0.5],'filled','XJitter','rand','XJitterWidth',0.5);
swarmchart([ones(numel(deltab_null),1)],abs(deltab_null),5,[0.5 0.5 0.5],'filled','XJitterWidth',0.5);
scatter(1,abs(deltab),20,[1 0 0],'filled')
set(gca,'xLim',[0.5 1.5],'xTick',[1], 'xTickLabel', {'Δbeta null'},'FontName', 'Arial');
set(gca,'yLim',[0 0.3],'yTick',[0 0.3], 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';

subplot(133); hold on;
% swarmchart([ones(numel(deltaR2_null),1)],deltaR2_null,5,[0.5 0.5 0.5],'filled','XJitter','rand','XJitterWidth',0.5);
swarmchart([ones(numel(bm_null),1)],abs(bm_null),5,[0.5 0.5 0.5],'filled','XJitterWidth',0.5);
scatter(1,abs(bm),20,[1 0 0],'filled')
set(gca,'xLim',[0.5 1.5],'xTick',[1], 'xTickLabel', {'ccsi beta null'},'FontName', 'Arial');
set(gca,'yLim',[0 0.6],'yTick',[0 0.6], 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';


%% Complex cortex: bb>0
% close all
% clc
mask = CBFsnr_mean_bil>4  &  bb_g>0 ; 

% plotting
figure('Color',[1 1 1],'Position',[0 0 300 300],'Units','pixels');
hold on;
% subfun_plot_scatter_density_isopotential(stru,func, 2 , 1 , 1 , 0 )
subfun_plot_scatter_regression(bb_g(mask),fc1_g(mask), [-0.12 0.12], [-8 8]);


% set
func = zscore(fc1_g(mask)); 
stru = zscore(bb_g(mask)); 
meta_cbf = zscore(cbf(mask)); 
null1 = zscore(cbf_null(mask,:));
meta_ccsi = zscore(ccsi(mask)); 
null2 = zscore(ccsi_null(mask,:));
fprintf('For bb gradient>0:\n\n');



% VIF
% 计算相关系数矩阵的逆，其对角线元素即为 VIF
% 或者使用统计工具箱自带函数
R = corrcoef([stru, meta_ccsi]);
vif_values = diag(inv(R))'; 
fprintf('VIF values: G_MPC = %.2f, CCSI = %.2f\n', vif_values(1), vif_values(2));




% ---------- Step 1: baseline model ----------
X1 = [ones(length(stru),1), stru];   % 截距 + Stru
[b1,~,~,~,stats1] = regress(func, X1);
R2_1 = stats1(1); % R^2
fprintf('Basic model:\nStep1: Func = %.3f + %.3f*Stru, R^2=%.3f, p=%.3g\n', b1(1), b1(2), R2_1, stats1(3));



% 第二步：可视化残差与 CCSI 的关系
% residuals_base = func - b1(1) - b1(2)*stru;
% figure;
% scatter(meta_ccsi, residuals_base, 'filled', 'MarkerFaceAlpha', 0.5);
% hold on;
% subfun_plot_scatter_regression(  meta_ccsi , residuals_base , [-3 3] , [-2 2]);
% xlabel('CCSI (Perfusion-Structure Coupling)');
% ylabel('Residuals from G_{MPC} model');
% title('Does CCSI explain the "unexplained" functional variance?');

% 第三步：量化 CCSI 对残差的解释力
% r_val = corr(meta_ccsi, residuals_base, 'type', 'Pearson');
% r_val_null = corr(null2, residuals_base, 'type', 'Pearson');
% p_val = mean(abs(r_val)<abs(r_val_null));
% fprintf('CCSI vs Residuals: r = %.3f, p = %.4f\n', r_val, p_val);
% r_val = corr(meta_ccsi, func, 'type', 'Pearson');
% r_val_null = corr(null2, func, 'type', 'Pearson');
% p_val = mean(abs(r_val)<abs(r_val_null));
% fprintf('CCSI vs func: r = %.3f, p = %.4f\n', r_val, p_val);




% ---------- Step 2: additional model (ccsi) ----------
X2 = [ones(length(stru),1), stru, meta_ccsi]; % 截距 + Stru + Meta
[b2,~,~,~,stats2] = regress(func, X2);
R2_2 = stats2(1); % R^2
fprintf('\nAdd CCSI:\nStep3: Func = %.3f + %.3f*Stru + %.3f*CCSI, R^2=%.3f, p=%.3g\n', ...
    b2(1), b2(2), b2(3), R2_2, stats2(3));

% spin-null test
nperm = size(null2,2);
deltaR2_null = zeros(nperm,1);
deltab_null = zeros(nperm,1);
bm_null = zeros(nperm,1);
for ip = 1:nperm
    X2_perm = [ones(length(stru),1), stru, null2(:,ip)];
    [bn,~,~,~,stats_perm] = regress(func, X2_perm);
    deltaR2_null(ip) = stats_perm(1) - R2_1;
    deltab_null(ip) = bn(2) - b1(2);
    bm_null(ip) = bn(3);
end

deltaR2 = R2_2 - R2_1; 
p_permr2 = mean(deltaR2_null >= deltaR2);
fprintf('ΔR² (Increased explanatory power) = %.3f,  ', deltaR2);
fprintf('ΔR² spin-null p-value = %.3g\n', p_permr2);

deltab = b2(2) - b1(2);
p_permb = mean(abs(deltab_null) >= abs(deltab));
fprintf('Δbeta (increased stru-func coupling) = %.3f,  ', deltab);
fprintf('Δbeta spin-null p-value = %.3g\n', p_permb);

allp(5) = p_permr2;
allp(6) = p_permb;

bm = b2(3);
p_permb = mean(abs(bm_null) >= abs(bm));
fprintf('ccsi beta = %.3f,  ', bm);
fprintf('ccsi beta spin-null p-value = %.3g\n', p_permb);

allp(9) = p_permb;



% plot
figure('Color',[1 1 1],'Position',[0 0 600 300],'Units','pixels');
subplot(141); hold on;
bar(1, R2_1, 'FaceColor', 'none', 'EdgeColor', [0, 200, 180]/255,'linewidth',1);
bar(3, R2_2, 'FaceColor', 'none', 'EdgeColor', [150, 103, 185]/255,'linewidth',1);
set(gca,'xLim',[0 4],'xTick',[1 3], 'xTickLabel', {'basic R2','add CCSI R2'},'FontName', 'Arial');
set(gca,'yLim',[0 0.2],'yTick',[0 0.2], 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';

subplot(142); hold on;
% swarmchart([ones(numel(deltaR2_null),1)],deltaR2_null,5,[0.5 0.5 0.5],'filled','XJitter','rand','XJitterWidth',0.5);
swarmchart([ones(numel(deltaR2_null),1)],deltaR2_null,5,[0.5 0.5 0.5],'filled','XJitterWidth',0.5);
scatter(1,deltaR2,20,[1 0 0],'filled')
set(gca,'xLim',[0.5 1.5],'xTick',[1], 'xTickLabel', {'ΔR2 null'},'FontName', 'Arial');
set(gca,'yLim',[0 0.4],'yTick',[0 0.4], 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';

subplot(143); hold on;
bar(1, b1(2), 'FaceColor', 'none', 'EdgeColor', [0, 200, 180]/255,'linewidth',1);
bar(3, b2(2), 'FaceColor', 'none', 'EdgeColor', [150, 103, 185]/255,'linewidth',1);
set(gca,'xLim',[0 4],'xTick',[1 3], 'xTickLabel', {'basic beta','add CCSI beta'},'FontName', 'Arial');
set(gca,'yLim',[0 0.3],'yTick',[0 0.3], 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';

subplot(144); hold on;
% swarmchart([ones(numel(deltab_null),1)],abs(deltab_null),5,[0.5 0.5 0.5],'filled','XJitter','rand','XJitterWidth',0.5);
swarmchart([ones(numel(deltab_null),1)],abs(deltab_null),5,[0.5 0.5 0.5],'filled','XJitterWidth',0.5);
scatter(1,deltab,20,[1 0 0],'filled')
set(gca,'xLim',[0.5 1.5],'xTick',[1], 'xTickLabel', {'Δbeta null'},'FontName', 'Arial');
set(gca,'yLim',[0 0.2],'yTick',[0 0.2], 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';







% plot 20260406
figure('Color',[1 1 1],'Position',[0 0 600 300],'Units','pixels');
subplot(131); hold on;
% swarmchart([ones(numel(deltaR2_null),1)],deltaR2_null,5,[0.5 0.5 0.5],'filled','XJitter','rand','XJitterWidth',0.5);
swarmchart([ones(numel(deltaR2_null),1)],deltaR2_null,5,[0.5 0.5 0.5],'filled','XJitterWidth',0.5);
scatter(1,deltaR2,20,[1 0 0],'filled')
set(gca,'xLim',[0.5 1.5],'xTick',[1], 'xTickLabel', {'ΔR2 null'},'FontName', 'Arial');
set(gca,'yLim',[0 0.4],'yTick',[0 0.4], 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';

subplot(132); hold on;
% swarmchart([ones(numel(deltaR2_null),1)],deltaR2_null,5,[0.5 0.5 0.5],'filled','XJitter','rand','XJitterWidth',0.5);
swarmchart([ones(numel(deltab_null),1)],abs(deltab_null),5,[0.5 0.5 0.5],'filled','XJitterWidth',0.5);
scatter(1,abs(deltab),20,[1 0 0],'filled')
set(gca,'xLim',[0.5 1.5],'xTick',[1], 'xTickLabel', {'Δbeta null'},'FontName', 'Arial');
set(gca,'yLim',[0 0.2],'yTick',[0 0.2], 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';

subplot(133); hold on;
% swarmchart([ones(numel(deltaR2_null),1)],deltaR2_null,5,[0.5 0.5 0.5],'filled','XJitter','rand','XJitterWidth',0.5);
swarmchart([ones(numel(bm_null),1)],abs(bm_null),5,[0.5 0.5 0.5],'filled','XJitterWidth',0.5);
scatter(1,abs(bm),20,[1 0 0],'filled')
set(gca,'xLim',[0.5 1.5],'xTick',[1], 'xTickLabel', {'ccsi beta null'},'FontName', 'Arial');
set(gca,'yLim',[0 0.6],'yTick',[0 0.6], 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';


%% FDR
FDR = fdr_BH(allp(1:9),0.05)



%% revision partial corr
% close all
% clc
mask = CBFsnr_mean_bil>4  &  bb_g>0 ; 


% set
func = zscore(fc1_g(mask)); 
stru = zscore(bb_g(mask)); 
meta_cbf = zscore(cbf(mask)); 
null1 = zscore(cbf_null(mask,:));
meta_ccsi = zscore(ccsi(mask)); 
null2 = zscore(ccsi_null(mask,:));
fprintf('For bb gradient>0:\n\n');



% plotting
figure('Color',[1 1 1],'Position',[0 0 300 300],'Units','pixels');
% hold on;
% % subfun_plot_scatter_density_isopotential(stru,func, 2 , 1 , 1 , 0 )
% subfun_plot_scatter_regression(stru, stru, [-2.5 2.5], [-2.5 2.5]);
subplot(131)
scatter(stru , func , 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
subfun_plot_scatter_regression(  stru , func , [-2.5 2.5] , [-2.5 2.5]);
xlabel('G_{MPC}');
ylabel('G_{func}');
title('G_{MPC} - G_{func}');



% VIF
% 计算相关系数矩阵的逆，其对角线元素即为 VIF
% 或者使用统计工具箱自带函数
R = corrcoef([stru, meta_ccsi]);
vif_values = diag(inv(R))'; 
fprintf('VIF values: G_MPC = %.2f, CCSI = %.2f\n', vif_values(1), vif_values(2));




% ---------- Step 1: baseline model ----------
X1 = [ones(length(stru),1), stru];   % 截距 + Stru
[b1,~,~,~,stats1] = regress(func, X1);
R2_1 = stats1(1); % R^2
fprintf('Basic model:\nStep1: Func = %.3f + %.3f*Stru, R^2=%.3f, p=%.3g\n', b1(1), b1(2), R2_1, stats1(3));


% 
% % 第二步：可视化残差与 CCSI 的关系
% residuals_base = func - b1(1) - b1(2)*stru;
% figure;
% scatter(meta_ccsi, residuals_base, 'filled', 'MarkerFaceAlpha', 0.5);
% hold on;
% subfun_plot_scatter_regression(  meta_ccsi , residuals_base , [-3 3] , [-2 2]);
% xlabel('CCSI (Perfusion-Structure Coupling)');
% ylabel('Residuals from G_{MPC} model');
% title('Does CCSI explain the "unexplained" functional variance?');
% 
% % 第三步：量化 CCSI 对残差的解释力
% r_val = corr(meta_ccsi, residuals_base, 'type', 'Pearson');
% r_val_null = corr(null2, residuals_base, 'type', 'Pearson');
% p_val = mean(abs(r_val)<abs(r_val_null));
% fprintf('CCSI vs Residuals: r = %.3f, p = %.4f\n', r_val, p_val);
% r_val = corr(meta_ccsi, func, 'type', 'Pearson');
% r_val_null = corr(null2, func, 'type', 'Pearson');
% p_val = mean(abs(r_val)<abs(r_val_null));
% fprintf('CCSI vs func: r = %.3f, p = %.4f\n', r_val, p_val);
% 
% residuals_base = func - b1(1) - b1(3)*meta_ccsi;
% figure;
% scatter(stru, residuals_base, 'filled', 'MarkerFaceAlpha', 0.5);
% hold on;
% subfun_plot_scatter_regression(  stru , residuals_base , [-3 3] , [-2 2]);
% xlabel('CCSI (Perfusion-Structure Coupling)');
% ylabel('Residuals from G_{MPC} model');
% title('Does CCSI explain the "unexplained" functional variance?');



% ---------- Step 2: additional model (ccsi) ----------
X2 = [ones(length(stru),1), stru, meta_ccsi]; % 截距 + Stru + Meta
[b2,~,~,~,stats2] = regress(func, X2);
R2_2 = stats2(1); % R^2
fprintf('\nAdd CCSI:\nStep3: Func = %.3f + %.3f*Stru + %.3f*CCSI, R^2=%.3f, p=%.3g\n', ...
    b2(1), b2(2), b2(3), R2_2, stats2(3));



% 第二步：可视化残差与 CCSI 的关系
residuals_base = func - b2(1) - b2(2)*stru;
% figure;
subplot(133)
scatter(meta_ccsi, residuals_base, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
subfun_plot_scatter_regression(  meta_ccsi , residuals_base , [-2.5 2.5] , [-2.5 2.5]);
xlabel('CCSI');
ylabel('G_{func} Residuals from G_{MPC}');
title('CCSI - Residuals(func-MPC)');

% 第三步：量化 CCSI 对残差的解释力
r_val = corr(meta_ccsi, residuals_base, 'type', 'Pearson');
r_val_null = corr(null2, residuals_base, 'type', 'Pearson');
p_val = mean(abs(r_val)<abs(r_val_null));
fprintf('CCSI vs Residuals: r = %.3f, p = %.4f\n', r_val, p_val);
r_val = corr(meta_ccsi, func, 'type', 'Pearson');
r_val_null = corr(null2, func, 'type', 'Pearson');
p_val = mean(abs(r_val)<abs(r_val_null));
fprintf('CCSI vs func: r = %.3f, p = %.4f\n', r_val, p_val);

residuals_base = func - b2(1) - b2(3)*meta_ccsi;
% figure;
subplot(132)
scatter(stru, residuals_base, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
subfun_plot_scatter_regression(  stru , residuals_base , [-2.5 2.5] , [-2.5 2.5]);
xlabel('G_{MPC}');
ylabel('G_{func} Residuals from CCSI');
title('G_{MPC} - Residuals(func-CCSI)');



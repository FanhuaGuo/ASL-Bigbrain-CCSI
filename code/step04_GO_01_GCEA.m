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
step01_CBFandCCSImaps_00_ComputeBasicResults;
step03_CellType_00_GeneInit;


%% mask
Sample_Mask = Sample_NumLable;
mask = CBFsnr_mean_uni>0;
GeneMask = ones(size(Sample_NumLable));
for i = 1:numel(Sample_Mask)
    if ~mask(Sample_Mask(i))
        GeneMask(i) = 0;
        Sample_Mask(i) = 0;
    end
end
Sample_Mask(Sample_Mask==0) = [];
GeneMask = GeneMask>0;
Sample_Mask_SNR0 = Sample_Mask;
GeneMask_SNR0 = GeneMask;

Sample_Mask = Sample_NumLable;
mask = CBFsnr_mean_uni>4;
GeneMask = ones(size(Sample_NumLable));
for i = 1:numel(Sample_Mask)
    if ~mask(Sample_Mask(i))
        GeneMask(i) = 0;
        Sample_Mask(i) = 0;
    end
end
Sample_Mask(Sample_Mask==0) = [];
GeneMask = GeneMask>0;
Sample_Mask_SNR4 = Sample_Mask;
GeneMask_SNR4 = GeneMask;


%% set 
CorrType = 'Pearson';  % Pearson Spearman

GeneNum_GO = [];
for i = 1:size(GO_categories,1)
    GeneNum_GO = [GeneNum_GO; numel(GO_categories.gene_label{i})];
end


%% ======================= Method1 median correlation's r-fisherZ ========================
%  =======================================================================================
%% common parameters
UniOrBil = '_uni';
BH_alpha = 0.05;
Strategy = 2; % 1:direct mean;  2:direct median
cate = GO_categories.gene_label;
nullType = '_variogram';  %  _variogram  _spin  _eigenstrapping


%% compute CBF CellType1-GSEA by median correlation's r-fisherZ
DataName = ['CBF_PC1_Val_L' num2str(num_layer)];
GeneMask = GeneMask_SNR0;
Sample_Mask = Sample_Mask_SNR0;
% nullType = '_variogram';  %  _variogram  _spin  _eigenstrapping


% read data
RealData_path = [Group_dir DataName UniOrBil '.1D'];
NullData_path = [Group_dir DataName UniOrBil nullType '.1D'];
GeneData = GeneExpression_Orig(:,GeneMask);
MapData = load(RealData_path);
MapData_null = load(NullData_path);
RepeatNum = size(MapData_null,2);


% compute GSEA
r = corr(MapData(Sample_Mask),GeneData','Type',CorrType, 'Rows', 'pairwise');
r_Fz = subfun_fisher_z(r);
rn = corr(MapData_null(Sample_Mask,:),GeneData','Type',CorrType, 'Rows', 'pairwise');
rn_Fz = subfun_fisher_z(rn);
cate_rZ = [];
cate_rnZ = [];
for i = 1:numel(cate)
    switch Strategy
        case 1
            cate_rZ = [cate_rZ mean(r_Fz(cate{i}))];
            cate_rnZ = [cate_rnZ mean(rn_Fz(:,cate{i}),2)];
        case 2
            cate_rZ = [cate_rZ median(r_Fz(cate{i}))];
            cate_rnZ = [cate_rnZ median(rn_Fz(:,cate{i}),2)];
    end
end
z = (cate_rZ-mean(cate_rnZ,1))./std(cate_rnZ);
spin_p = 2 * (1 - normcdf(abs(z)));
% spin_p = (0.01+sum(abs(cate_rZ) < abs(cate_rnZ))) / (0.01+RepeatNum);
BHp = fdr_BH(spin_p,BH_alpha)';


% disposal result
GCEAresults = table(GO_categories.go_id, GO_categories.go_func, cate_rZ', spin_p', BHp, GeneNum_GO,...
                    'VariableNames',{'go_id', 'go_func', 'GCEAscore', 'spin_p', 'FDR', 'GeneNum'});

GCEA_CBFscore_1 = GCEAresults;


z = (r-mean(rn,1))./std(rn);
spin_p = 2 * (1 - normcdf(abs(z)));
% spin_p = (1+sum(abs(r) < abs(rn))) / (1+RepeatNum);
BHp = fdr_BH(spin_p,BH_alpha)';
GeneResults = table(GeneID, r', spin_p', BHp,...
                    'VariableNames',{'gene_id', 'r', 'spin_p', 'FDR'});
GeneCorr_CBFscore_1 = GeneResults;



%% compute CCSI CellType1-GSEA by median correlation's r-fisherZ
DataName = ['CBF_Mean_PScs_L' num2str(num_layer)];
GeneMask = GeneMask_SNR4;
Sample_Mask = Sample_Mask_SNR4;
% nullType = '_eigenstrapping';  %  _variogram  _spin  _eigenstrapping


% read data
RealData_path = [Group_dir DataName UniOrBil '.1D'];
NullData_path = [Group_dir DataName UniOrBil nullType '.1D'];
GeneData = GeneExpression_Orig(:,GeneMask);
MapData = load(RealData_path);
MapData_null = load(NullData_path);
RepeatNum = size(MapData_null,2);


% compute GSEA
r = corr(MapData(Sample_Mask),GeneData','Type',CorrType, 'Rows', 'pairwise');
r_Fz = subfun_fisher_z(r);
rn = corr(MapData_null(Sample_Mask,:),GeneData','Type',CorrType, 'Rows', 'pairwise');
rn_Fz = subfun_fisher_z(rn);
cate_rZ = [];
cate_rnZ = [];
for i = 1:numel(cate)
    switch Strategy
        case 1
            cate_rZ = [cate_rZ mean(r_Fz(cate{i}))];
            cate_rnZ = [cate_rnZ mean(rn_Fz(:,cate{i}),2)];
        case 2
            cate_rZ = [cate_rZ median(r_Fz(cate{i}))];
            cate_rnZ = [cate_rnZ median(rn_Fz(:,cate{i}),2)];
    end
end
z = (cate_rZ-mean(cate_rnZ,1))./std(cate_rnZ);
spin_p = 2 * (1 - normcdf(abs(z)));
% spin_p = (0.01+sum(abs(cate_rZ) < abs(cate_rnZ))) / (0.01+RepeatNum);
BHp = fdr_BH(spin_p,BH_alpha)';


% disposal result
GCEAresults = table(GO_categories.go_id, GO_categories.go_func, cate_rZ', spin_p', BHp, GeneNum_GO,...
                    'VariableNames',{'go_id', 'go_func', 'GCEAscore', 'spin_p', 'FDR', 'GeneNum'});

GCEA_CCSI_1 = GCEAresults;



z = (r-mean(rn,1))./std(rn);
spin_p = 2 * (1 - normcdf(abs(z)));
% spin_p = (1+sum(abs(r) < abs(rn))) / (1+RepeatNum);
BHp = fdr_BH(spin_p,BH_alpha)';
GeneResults = table(GeneID, r', spin_p', BHp,...
                    'VariableNames',{'gene_id', 'r', 'spin_p', 'FDR'});
GeneCorr_CCSI_1 = GeneResults;


%% compute CCSI-control CellType1-GSEA by median correlation's r-fisherZ
% 'T1vCBF_Mean_PScs_L', 'T1vBB_Mean_PScs_L', 'T1wCBF_Mean_PScs_L', 'T1wBB_Mean_PScs_L'
DataName = ['T1vBB_Mean_PScs_L' num2str(num_layer)];
GeneMask = GeneMask_SNR4;
Sample_Mask = Sample_Mask_SNR4;
% nullType = '_variogram';  %  _variogram  _spin  _eigenstrapping


% read data
RealData_path = [Group_dir DataName UniOrBil '.1D'];
NullData_path = [Group_dir DataName UniOrBil nullType '.1D'];
GeneData = GeneExpression_Orig(:,GeneMask);
MapData = load(RealData_path);
MapData_null = load(NullData_path);
RepeatNum = size(MapData_null,2);


% compute GSEA
r = corr(MapData(Sample_Mask),GeneData','Type',CorrType, 'Rows', 'pairwise');
r_Fz = subfun_fisher_z(r);
rn = corr(MapData_null(Sample_Mask,:),GeneData','Type',CorrType, 'Rows', 'pairwise');
rn_Fz = subfun_fisher_z(rn);
cate_rZ = [];
cate_rnZ = [];
for i = 1:numel(cate)
    switch Strategy
        case 1
            cate_rZ = [cate_rZ mean(r_Fz(cate{i}))];
            cate_rnZ = [cate_rnZ mean(rn_Fz(:,cate{i}),2)];
        case 2
            cate_rZ = [cate_rZ median(r_Fz(cate{i}))];
            cate_rnZ = [cate_rnZ median(rn_Fz(:,cate{i}),2)];
    end
end
z = (cate_rZ-mean(cate_rnZ,1))./std(cate_rnZ);
spin_p = 2 * (1 - normcdf(abs(z)));
% spin_p = (0.01+sum(abs(cate_rZ) < abs(cate_rnZ))) / (0.01+RepeatNum);
BHp = fdr_BH(spin_p,BH_alpha)';


% disposal result
GCEAresults = table(GO_categories.go_id, GO_categories.go_func, cate_rZ', spin_p', BHp, GeneNum_GO,...
                    'VariableNames',{'go_id', 'go_func', 'GCEAscore', 'spin_p', 'FDR', 'GeneNum'});

GCEA_CCSIcon_1 = GCEAresults;


z = (r-mean(rn,1))./std(rn);
spin_p = 2 * (1 - normcdf(abs(z)));
% spin_p = (0.01+sum(abs(r) < abs(rn))) / (0.01+RepeatNum);
BHp = fdr_BH(spin_p,BH_alpha)';
GeneResults = table(GeneID, r', spin_p', BHp,...
                    'VariableNames',{'gene_id', 'r', 'spin_p', 'FDR'});
GeneCorr_CCSIcon_1 = GeneResults;




%% plotting volcano plot - CCSI
T = GeneCorr_CCSI_1;  % GeneCorr_CCSI_1 GeneCorr_CBFscore_1

% 如果 gene_id 被读成了带引号的字符串，也没关系
gene = string(T.gene_id);

% 2) 定义火山图坐标
% 横轴：r（如果你有 log2FC，建议换成 log2FC）
x = T.r;

% 纵轴：-log10(FDR)
% y = T.FDR;
y = -log10(T.FDR);

% 3) 阈值设置
xCutoff = 0;      % 你可以改成 1、0.5 等
fdrCutoff = 0.05;

% mask = abs(x) > xCutoff;
% y(mask) = fdr_BH(T.spin_p(mask),BH_alpha)';
% y = -log10(y);

% 4) 分类
isUp   = (x >=  xCutoff) & (y > -log10(fdrCutoff));
isDown = (x <= -xCutoff) & (y > -log10(fdrCutoff));
isNS   = ~(isUp | isDown);

% 5) 画图
figure('Color',[1 1 1],'Position',[0 0 300 400],'Units','pixels');
hold on;

scatter(x(isNS),   y(isNS),   18, [0.75 0.75 0.75], 'filled');
scatter(x(isUp),   y(isUp),   22, [0.85 0.20 0.20], 'filled');
scatter(x(isDown), y(isDown), 22, [0.20 0.35 0.85], 'filled');

% 阈值线
xline(xCutoff,  '--k', 'LineWidth', 1);
xline(-xCutoff, '--k', 'LineWidth', 1);
yline(-log10(fdrCutoff), '--k', 'LineWidth', 1);

% 6) 标注最显著的几个基因
% 这里按 FDR 从小到大取前 10 个
% [~, idx] = sort(T.FDR, 'ascend');
% topN = min(10, height(T));
% idxTop = idx(1:topN);
% 
% text(x(idxTop), y(idxTop), " " + gene(idxTop), ...
%     'FontSize', 9, 'Color', 'k', 'Interpreter', 'none');

% 7) 美化
set(gca, 'xlim', [-1 1], 'xtick', [-1:0.5:1]);
set(gca, 'ylim', [0 4], 'ytick', [0:4], 'yTickLabel', {'1','0.1','0.01','0.001','0.0001'});
% xlabel('r');
% ylabel('-log_{10}(FDR)');
% title('Volcano Plot');
grid on;
box on;

% legend({'Not significant','Up','Down'}, 'Location', 'best');
set(gca, 'FontSize', 12);







%% plotting GO - CCSI
ifout = 0;
PlotData = GCEA_CCSI_1;


% simple process
PlotData = sortrows(PlotData,5);
PlotData = PlotData(PlotData.FDR<0.05,:);
if ifout
    tmp = PlotData(:,1);
    writetable(tmp, [Group_dir 'GO_sigCCSI_goid.csv']);
end
tmp = readtable([Group_dir 'GO_sigCCSI_goid_wCate.csv']);
PlotData = [PlotData tmp(:,3)];

x = PlotData.GCEAscore;
y = -log10(PlotData.FDR);
tmpcate = PlotData.category;
isDM = zeros(numel(tmpcate),1);
isIM = zeros(numel(tmpcate),1);
isNM = zeros(numel(tmpcate),1);
for i = 1:numel(tmpcate)
    switch tmpcate{i}
        case 'direct_metabolism'
            isDM(i) = 1;
        case 'indirect_metabolism'
            isIM(i) = 1;
        case 'non_metabolism'
            isNM(i) = 1;
    end
end
isDM = isDM>0;
isIM = isIM>0;
isNM = isNM>0;


% 5) 画图
figure('Color',[1 1 1],'Position',[0 0 600 400],'Units','pixels');
hold on;

scatter(x(isNM), y(isNM), 200, [0.75 0.75 0.75], 'filled', 'MarkerEdgeColor', 'k');
scatter(x(isIM), y(isIM), 200, [0.20 0.35 0.85], 'filled', 'MarkerEdgeColor', 'k');
scatter(x(isDM), y(isDM), 200, [0.85 0.20 0.20], 'filled', 'MarkerEdgeColor', 'k');

% 阈值线
% xline(xCutoff,  '--k', 'LineWidth', 1);
% xline(-xCutoff, '--k', 'LineWidth', 1);
% yline(-log10(fdrCutoff), '--k', 'LineWidth', 1);

% 6) 标注最显著的几个基因
% 这里按 FDR 从小到大取前 10 个
% [~, idx] = sort(T.FDR, 'ascend');
% topN = min(10, height(T));
% idxTop = idx(1:topN);
% 
% text(x(idxTop), y(idxTop), " " + gene(idxTop), ...
%     'FontSize', 9, 'Color', 'k', 'Interpreter', 'none');

% 7) 美化
set(gca, 'xlim', [-0.6 0.6], 'xtick', [-1:0.5:1]);
set(gca, 'ylim', [1 4], 'ytick', [0:1 -log10(0.05) 2:4], 'yTickLabel', {'1','0.1','0.05','0.01','0.001','0.0001'});
% xlabel('r');
% ylabel('-log_{10}(FDR)');
grid off;
box on;

% legend({'Not significant','Up','Down'}, 'Location', 'best');
set(gca, 'FontSize', 12);




%%






%% plotting CBF old
threshold = -log10([0.05 0.01]);
Plotdata = GCEA_CBFscore_1;
Plotdata = Plotdata(Plotdata.FDR <= 0.05,:);
Plotdata = sortrows(Plotdata,'GCEAscore','ascend');
Plotdata(Plotdata.GCEAscore<0,:) = [];


% classification color
classcolor = [255 0   0;      % development and differentiation
              150 150 150;    % others
              102 45  145;    % metabolism
              241 90  41;     % vascular
              28  117 188;    % stress response and immune signaling
              71  0   184;    % signal transduction / signaling regulation
              255 255 255;    % basic process
              ]./255;
Plotclass = [5 1 1 1 1 ...
             1 1 1];


% bar color
Plotscore = Plotdata.GCEAscore;
Plotp = -log10(Plotdata.FDR);
Plotp(Plotp>threshold(2)) = threshold(2);
Plotp(Plotp<threshold(1)) = threshold(1);
Plotp = 1+round(1000 * (Plotp - threshold(1)) / (threshold(2) - threshold(1)));


% color gradient
colorGra = [ones(1,1001); 0.9:-0.9/1000:0; 0.9:-0.9/1000:0]';


% plotting
figure('Color',[1 1 1],'Position',[0 0 600 300],'Units','pixels');
hold on;
labels = Plotdata.go_func;
for i = 1:length(Plotscore)
    plot([0 Plotscore(i)], [i i], 'Color', colorGra(Plotp(i),:), 'linewidth', 2);
    plot(Plotscore(i), i, 'Marker', 'o', 'MarkerSize', 15, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', classcolor(Plotclass(i),:), 'linewidth', 0.5);
end
set(gca,'xLim',[0 0.45],'xTick',[0 0.4],'FontName', 'Arial');
set(gca,'yLim',[0 numel(Plotscore)+1],'yTick',1:numel(Plotscore),'yTickLabel',labels, 'FontName', 'Arial');
set(gca, 'Box', 'on');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';




%% plotting CCSI
threshold = -log10([0.05 0.001]);
Plotdata = GCEA_CCSI_1;
Plotdata = Plotdata(Plotdata.FDR <= 0.05,:);
Plotdata = sortrows(Plotdata,'GCEAscore','ascend');
Plotdata(Plotdata.GCEAscore<0,:) = [];


% classification color
classcolor = [255 0   0;      % development and differentiation
              150 150 150;    % others
              102 45  145;    % metabolism
              201 120  90;     % vascular
              28  117 188;    % stress response and immune signaling
              71  0   184;    % signal transduction / signaling regulation
              255 255 255;    % basic process
              ]./255;
Plotclass = [5 5 6 2 6 ...
             1 6 1 2 2 ...
             4 7 7 3 2 ...
             3 1 3 6 1 ...
             1 1 1 4 3 ...
             3 5 3 4 1 ...
             1 3]; 
% note 3,21,28,29,30 may link to NVC
% note 10,16,17,18,22,27 may link to MRC

% bar color
Plotscore = Plotdata.GCEAscore;
Plotp = -log10(Plotdata.FDR);
Plotp(Plotp>threshold(2)) = threshold(2);
Plotp(Plotp<threshold(1)) = threshold(1);
Plotp = 1+round(1000 * (Plotp - threshold(1)) / (threshold(2) - threshold(1)));


% color gradient
colorGra = [ones(1,1001); 0.9:-0.9/1000:0; 0.9:-0.9/1000:0]';


% plotting
figure('Color',[1 1 1],'Position',[0 0 600 700],'Units','pixels');
hold on;
labels = Plotdata.go_func;
for i = 1:length(Plotscore)
    plot([0 Plotscore(i)], [i i], 'Color', colorGra(Plotp(i),:), 'linewidth', 2);
    plot(Plotscore(i), i, 'Marker', 'o', 'MarkerSize', 15, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', classcolor(Plotclass(i),:), 'linewidth', 0.5);
end
set(gca,'xLim',[0 0.45],'xTick',[0 0.4],'FontName', 'Arial');
set(gca,'yLim',[0 numel(Plotscore)+1],'yTick',1:numel(Plotscore),'yTickLabel',labels, 'FontName', 'Arial');
set(gca, 'Box', 'on');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';








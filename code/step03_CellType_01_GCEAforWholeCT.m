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


%% change the plotting order
newOrder = [2 25 3:21 1 22 23 24];
CT_categories1 = CT_categories1(newOrder,:);


%% set 
CorrType = 'Pearson';  % Pearson Spearman


%% ======================= Method1 median correlation's r-fisherZ ========================
%  =======================================================================================
%% common parameters
UniOrBil = '_uni';
BH_alpha = 0.05;
Strategy = 2; % 1:direct mean;  2:direct median
cate = CT_categories1.gene_label;
cateLabel = CT_categories1.Cell;
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
spin_p = (1+sum(abs(cate_rZ) < abs(cate_rnZ))) / (1+RepeatNum);
BHp = fdr_BH(spin_p,BH_alpha)';


% disposal result
GCEAresults = table(cateLabel, cate_rZ', spin_p', BHp,...
                    'VariableNames',{'Cell', 'GCEAscore', 'spin_p', 'FDR'});

GCEA_CBFscore_1 = GCEAresults;
                
      


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
spin_p = (1+sum(abs(cate_rZ) < abs(cate_rnZ))) / (1+RepeatNum);
BHp = fdr_BH(spin_p,BH_alpha)';


% disposal result
GCEAresults = table(cateLabel, cate_rZ', spin_p', BHp,...
                    'VariableNames',{'Cell', 'GCEAscore', 'spin_p', 'FDR'});

GCEA_CCSI_1 = GCEAresults;

      


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
spin_p = (1+sum(abs(cate_rZ) < abs(cate_rnZ))) / (1+RepeatNum);
BHp = fdr_BH(spin_p,BH_alpha)';


% disposal result
GCEAresults = table(cateLabel, cate_rZ', spin_p', BHp,...
                    'VariableNames',{'Cell', 'GCEAscore', 'spin_p', 'FDR'});

GCEA_CCSIcon_1 = GCEAresults;




%% plotting CBF
% manual parameters
DataPlot = GCEA_CBFscore_1;


% auto parameters
cellLabel = DataPlot.Cell;
PlotScore = DataPlot.GCEAscore;
PlotFDR = DataPlot.FDR;
Plotp = DataPlot.spin_p;


% plotting
figure('Color',[1 1 1],'Position',[0 0 500 800],'Units','pixels');
hold on;
y = numel(cellLabel):-1:1;
for i = 1:numel(PlotScore)
    if PlotFDR(i) <= 0.05
        plot([0 PlotScore(i)], [y(i) y(i)], 'r--','LineWidth',1);
        plot(PlotScore(i), y(i), 'ro','MarkerFaceColor','r','MarkerSize',6);
    elseif Plotp(i) <= 0.05
        plot([0 PlotScore(i)], [y(i) y(i)], 'k--','LineWidth',1);
        plot(PlotScore(i), y(i), 'ko','MarkerFaceColor','k','MarkerSize',6);
    else
        plot([0 PlotScore(i)], [y(i) y(i)], 'k--','LineWidth',1);
        plot(PlotScore(i), y(i), 'ko','MarkerFaceColor','k','MarkerSize',6);
    end
end
plot([0 0],[0 numel(PlotScore)+1],'-','color',[0.7 0.7 0.7],'LineWidth',1);
set(gca,'xLim',[-0.4 0.4],'xTick',[-0.3 0 0.3],'FontName', 'Arial');
set(gca,'yLim',[0.5 numel(PlotScore)+0.5],'yTick',[1:numel(PlotScore)],'yTickLabel',cellLabel(y), 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';



%% plotting CCSI
% manual parameters
DataPlot = GCEA_CCSI_1;
DataPlotcon = GCEA_CCSIcon_1;


% auto parameters
cellLabel = DataPlot.Cell;
PlotScore = DataPlot.GCEAscore;
PlotFDR = DataPlot.FDR;
Plotp = DataPlot.spin_p;


% plotting v2
% figure('Color',[1 1 1],'Position',[0 0 500 800],'Units','pixels');
% hold on;
% y = numel(cellLabel):-1:1;
% for i = 1:numel(PlotScore)
%     if PlotFDR(i) <= 0.05
%         plot([0 PlotScore(i)], [y(i) y(i)], 'r--','LineWidth',1);
%         plot(PlotScore(i), y(i), 'ro','MarkerFaceColor','r','MarkerSize',6);
%     elseif Plotp(i) <= 0.05
%         plot([0 PlotScore(i)], [y(i) y(i)], 'k--','LineWidth',1);
%         plot(PlotScore(i), y(i), 'ko','MarkerFaceColor','k','MarkerSize',6);
%     else
%         plot([0 PlotScore(i)], [y(i) y(i)], 'k--','LineWidth',1);
%         plot(PlotScore(i), y(i), 'ko','MarkerFaceColor','k','MarkerSize',6);
%     end
% end
% plot([0 0],[0 numel(PlotScore)+1],'-','color',[0.7 0.7 0.7],'LineWidth',1);
% set(gca,'xLim',[-0.4 0.4],'xTick',[-0.3 0 0.3],'FontName', 'Arial');
% set(gca,'yLim',[0.5 numel(PlotScore)+0.5],'yTick',[1:numel(PlotScore)],'yTickLabel',cellLabel(y), 'FontName', 'Arial');
% set(gca, 'Box', 'off');
% ax = gca;
% ax.YColor = 'k'; 
% ax.XColor = 'k'; % 保持刻度数字
% ax.TickDir = 'out';


% plotting v3
figure('Color',[1 1 1],'Position',[0 0 500 800],'Units','pixels');
hold on;
shift = 0.15;
y = numel(cellLabel):-1:1;
PlotFDR = DataPlotcon.FDR;
PlotScore = DataPlotcon.GCEAscore;
for i = 1:numel(PlotScore)
    if PlotFDR(i) <= 0.05
        plot([0 PlotScore(i)], [y(i)-shift y(i)-shift], 'r--','LineWidth',1);
        plot(PlotScore(i), y(i)-shift, 'rv','MarkerFaceColor','r','MarkerSize',6);
    else
        plot([0 PlotScore(i)], [y(i)-shift y(i)-shift], 'k--','LineWidth',1);
        plot(PlotScore(i), y(i)-shift, 'kv','MarkerFaceColor','k','MarkerSize',6);
    end
end
PlotFDR = DataPlot.FDR;
PlotScore = DataPlot.GCEAscore;
for i = 1:numel(PlotScore)
    if PlotFDR(i) <= 0.05
        plot([0 PlotScore(i)], [y(i)+shift y(i)+shift], 'r--','LineWidth',1);
        plot(PlotScore(i), y(i)+shift, 'r^','MarkerFaceColor','r','MarkerSize',6);
    else
        plot([0 PlotScore(i)], [y(i)+shift y(i)+shift], 'k--','LineWidth',1);
        plot(PlotScore(i), y(i)+shift, 'k^','MarkerFaceColor','k','MarkerSize',6);
    end
end
plot([0 0],[0 numel(PlotScore)+1],'-','color',[0.7 0.7 0.7],'LineWidth',1);
set(gca,'xLim',[-0.4 0.4],'xTick',[-0.3 0 0.3],'FontName', 'Arial');
set(gca,'yLim',[0.5 numel(PlotScore)+0.5],'yTick',[1:numel(PlotScore)],'yTickLabel',cellLabel(y), 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';



%% ======================= Method2 median Expression maps ================================
%  =======================================================================================
%% compute cell distribution by median GeneExp
% common parameters
cate = CT_categories1.gene_label;
cateLabel = CT_categories1.Cell;
UniOrBil = '_uni';
BH_alpha = 0.05;

% compute
Cell_dis_m2 = [];
GeneData = GeneExpression_Orig;
for i = 1:numel(cate)
    Cell_dis_m2 = [Cell_dis_m2; median(GeneExpression_Orig(cate{i},:),1)];
end
Cell_dis_m2 = Cell_dis_m2';


%% compute CBF CellType1-GSEA by median GeneExp as Cell-distribution
DataName = ['CBF_PC1_Val_L' num2str(num_layer)];
GeneMask = GeneMask_SNR0;
Sample_Mask = Sample_Mask_SNR0;
% nullType = '_variogram';  %  _variogram  _spin  _eigenstrapping


% read data
RealData_path = [Group_dir DataName UniOrBil '.1D'];
NullData_path = [Group_dir DataName UniOrBil nullType '.1D'];
MapData = load(RealData_path);
MapData_null = load(NullData_path);
RepeatNum = size(MapData_null,2);


% compute GSEA
r = corr(MapData(Sample_Mask),Cell_dis_m2(GeneMask,:),'Type',CorrType, 'Rows', 'pairwise');
r_Fz = subfun_fisher_z(r);
rn = corr(MapData_null(Sample_Mask,:),Cell_dis_m2(GeneMask,:),'Type',CorrType, 'Rows', 'pairwise');
rn_Fz = subfun_fisher_z(rn);
spin_p = (1+sum(abs(r_Fz) < abs(rn_Fz)))/(1+RepeatNum);
BHp = fdr_BH(spin_p,BH_alpha)';


% disposal result
GCEAresults = table(cateLabel, r_Fz', spin_p', BHp,...
                    'VariableNames',{'Cell', 'GCEAscore', 'spin_p', 'FDR'});

GCEA_CBFscore_2 = GCEAresults;
                
      


%% compute CCSI CellType1-GSEA by median GeneExp as Cell-distribution
DataName = ['CBF_Mean_PScs_L' num2str(num_layer)];
GeneMask = GeneMask_SNR4;
Sample_Mask = Sample_Mask_SNR4;
% nullType = '_eigenstrapping';  %  _variogram  _spin  _eigenstrapping


% read data
RealData_path = [Group_dir DataName UniOrBil '.1D'];
NullData_path = [Group_dir DataName UniOrBil nullType '.1D'];
MapData = load(RealData_path);
MapData_null = load(NullData_path);
RepeatNum = size(MapData_null,2);


% compute GSEA
r = corr(MapData(Sample_Mask),Cell_dis_m2(GeneMask,:),'Type',CorrType, 'Rows', 'pairwise');
r_Fz = subfun_fisher_z(r);
rn = corr(MapData_null(Sample_Mask,:),Cell_dis_m2(GeneMask,:),'Type',CorrType, 'Rows', 'pairwise');
rn_Fz = subfun_fisher_z(rn);
spin_p = (1+sum(abs(r_Fz) < abs(rn_Fz)))/(1+RepeatNum);
BHp = fdr_BH(spin_p,BH_alpha)';


% disposal result
GCEAresults = table(cateLabel, r_Fz', spin_p', BHp,...
                    'VariableNames',{'Cell', 'GCEAscore', 'spin_p', 'FDR'});

GCEA_CCSI_2 = GCEAresults;




%% compute CCSI-control CellType1-GSEA by median correlation's r-fisherZ
% 'T1vCBF_Mean_PScs_L', 'T1vBB_Mean_PScs_L', 'T1wCBF_Mean_PScs_L', 'T1wBB_Mean_PScs_L'
DataName = ['T1vBB_Mean_PScs_L' num2str(num_layer)];
GeneMask = GeneMask_SNR4;
Sample_Mask = Sample_Mask_SNR4;
% nullType = '_eigenstrapping';  %  _variogram  _spin  _eigenstrapping


% read data
RealData_path = [Group_dir DataName UniOrBil '.1D'];
NullData_path = [Group_dir DataName UniOrBil nullType '.1D'];
MapData = load(RealData_path);
MapData_null = load(NullData_path);
RepeatNum = size(MapData_null,2);


% compute GSEA
r = corr(MapData(Sample_Mask),Cell_dis_m2(GeneMask,:),'Type',CorrType, 'Rows', 'pairwise');
r_Fz = subfun_fisher_z(r);
rn = corr(MapData_null(Sample_Mask,:),Cell_dis_m2(GeneMask,:),'Type',CorrType, 'Rows', 'pairwise');
rn_Fz = subfun_fisher_z(rn);
spin_p = (1+sum(abs(r_Fz) < abs(rn_Fz)))/(1+RepeatNum);
BHp = fdr_BH(spin_p,BH_alpha)';


% disposal result
GCEAresults = table(cateLabel, r_Fz', spin_p', BHp,...
                    'VariableNames',{'Cell', 'GCEAscore', 'spin_p', 'FDR'});

GCEA_CCSIcon_2 = GCEAresults;



%% write to SourceData Fig5cd EDFig6cd
DataName = ['CBF_PC1_Val_L' num2str(num_layer)];
outSD = load([Group_dir DataName UniOrBil '.1D']);
DataName = ['CBF_Mean_PScs_L' num2str(num_layer)];
outSD = [outSD load([Group_dir DataName UniOrBil '.1D'])];
DataName = ['T1vBB_Mean_PScs_L' num2str(num_layer)];
outSD = [outSD load([Group_dir DataName UniOrBil '.1D'])];
outSD = outSD(Sample_Mask,:);
outSD = [outSD Cell_dis_m2(GeneMask,[1 25 18 19])];
outSD2 = atlas_info.region(Sample_Mask);

%% update FDR
p = [GCEA_CBFscore_2.spin_p; GCEA_CCSI_2.spin_p; GCEA_CCSIcon_2.spin_p];
FDR = fdr_BH(p,0.5);
GCEA_CBFscore_2.FDRtog = FDR(1:end/3)';
GCEA_CCSI_2.FDRtog = FDR(end/3+1:2*end/3)';
GCEA_CCSIcon_2.FDRtog = FDR(2*end/3+1:end)';



%% plotting CBF
% manual parameters
DataPlot = GCEA_CBFscore_2;


% auto parameters
cellLabel = DataPlot.Cell;
PlotScore = DataPlot.GCEAscore;
PlotFDR = DataPlot.FDR;
Plotp = DataPlot.spin_p;


% plotting
figure('Color',[1 1 1],'Position',[0 0 500 800],'Units','pixels');
hold on;
y = numel(cellLabel):-1:1;
for i = 1:numel(PlotScore)
    if PlotFDR(i) <= 0.05
        plot([0 PlotScore(i)], [y(i) y(i)], 'r--','LineWidth',1);
        plot(PlotScore(i), y(i), 'ro','MarkerFaceColor','r','MarkerSize',6);
    elseif Plotp(i) <= 0.05
        plot([0 PlotScore(i)], [y(i) y(i)], 'k--','LineWidth',1);
        plot(PlotScore(i), y(i), 'ko','MarkerFaceColor','k','MarkerSize',6);
    else
        plot([0 PlotScore(i)], [y(i) y(i)], 'k--','LineWidth',1);
        plot(PlotScore(i), y(i), 'ko','MarkerFaceColor','k','MarkerSize',6);
    end
end
plot([0 0],[0 numel(PlotScore)+1],'-','color',[0.7 0.7 0.7],'LineWidth',1);
set(gca,'xLim',[-0.8 0.8],'xTick',[-0.6 0 0.6],'FontName', 'Arial');
set(gca,'yLim',[0.5 numel(PlotScore)+0.5],'yTick',[1:numel(PlotScore)],'yTickLabel',cellLabel(y), 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';



%% plotting CCSI v3
% manual parameters
DataPlot = GCEA_CCSI_2;
DataPlotcon = GCEA_CCSIcon_2;


% auto parameters
cellLabel = DataPlot.Cell;
PlotScore = DataPlot.GCEAscore;
PlotFDR = DataPlot.FDR;
Plotp = DataPlot.spin_p;


% plotting v2
% figure('Color',[1 1 1],'Position',[0 0 500 800],'Units','pixels');
% hold on;
% y = numel(cellLabel):-1:1;
% for i = 1:numel(PlotScore)
%     if PlotFDR(i) <= 0.05
%         plot([0 PlotScore(i)], [y(i) y(i)], 'r--','LineWidth',1);
%         plot(PlotScore(i), y(i), 'ro','MarkerFaceColor','r','MarkerSize',6);
%     elseif Plotp(i) <= 0.05
%         plot([0 PlotScore(i)], [y(i) y(i)], 'k--','LineWidth',1);
%         plot(PlotScore(i), y(i), 'ko','MarkerFaceColor','k','MarkerSize',6);
%     else
%         plot([0 PlotScore(i)], [y(i) y(i)], 'k--','LineWidth',1);
%         plot(PlotScore(i), y(i), 'ko','MarkerFaceColor','k','MarkerSize',6);
%     end
% end
% plot([0 0],[0 numel(PlotScore)+1],'-','color',[0.7 0.7 0.7],'LineWidth',1);
% set(gca,'xLim',[-0.4 0.4],'xTick',[-0.3 0 0.3],'FontName', 'Arial');
% set(gca,'yLim',[0.5 numel(PlotScore)+0.5],'yTick',[1:numel(PlotScore)],'yTickLabel',cellLabel(y), 'FontName', 'Arial');
% set(gca, 'Box', 'off');
% ax = gca;
% ax.YColor = 'k'; 
% ax.XColor = 'k'; % 保持刻度数字
% ax.TickDir = 'out';


% plotting v3
figure('Color',[1 1 1],'Position',[0 0 500 800],'Units','pixels');
hold on;
shift = 0.15;
y = numel(cellLabel):-1:1;
PlotFDR = DataPlotcon.FDR;
PlotScore = DataPlotcon.GCEAscore;
for i = 1:numel(PlotScore)
    if PlotFDR(i) <= 0.05
        plot([0 PlotScore(i)], [y(i)-shift y(i)-shift], 'r--','LineWidth',1);
        plot(PlotScore(i), y(i)-shift, 'rv','MarkerFaceColor','r','MarkerSize',6);
    else
        plot([0 PlotScore(i)], [y(i)-shift y(i)-shift], 'k--','LineWidth',1);
        plot(PlotScore(i), y(i)-shift, 'kv','MarkerFaceColor','k','MarkerSize',6);
    end
end
PlotFDR = DataPlot.FDR;
PlotScore = DataPlot.GCEAscore;
for i = 1:numel(PlotScore)
    if PlotFDR(i) <= 0.05
        plot([0 PlotScore(i)], [y(i)+shift y(i)+shift], 'r--','LineWidth',1);
        plot(PlotScore(i), y(i)+shift, 'r^','MarkerFaceColor','r','MarkerSize',6);
    else
        plot([0 PlotScore(i)], [y(i)+shift y(i)+shift], 'k--','LineWidth',1);
        plot(PlotScore(i), y(i)+shift, 'k^','MarkerFaceColor','k','MarkerSize',6);
    end
end
plot([0 0],[0 numel(PlotScore)+1],'-','color',[0.7 0.7 0.7],'LineWidth',1);
set(gca,'xLim',[-0.8 0.8],'xTick',[-0.6 0 0.6],'FontName', 'Arial');
set(gca,'yLim',[0.5 numel(PlotScore)+0.5],'yTick',[1:numel(PlotScore)],'yTickLabel',cellLabel(y), 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';



%% plotting CCSI(old)
% % manual parameters
% DataPlot = GCEA_CCSI_2;
% 
% 
% % auto parameters
% cellLabel = DataPlot.Cell;
% PlotScore = DataPlot.GCEAscore;
% PlotFDR = DataPlot.FDR;
% Plotp = DataPlot.spin_p;
% 
% 
% % plotting
% figure('Color',[1 1 1],'Position',[0 0 500 800],'Units','pixels');
% hold on;
% y = numel(cellLabel):-1:1;
% for i = 1:numel(PlotScore)
%     if PlotFDR(i) <= 0.05
%         plot([0 PlotScore(i)], [y(i) y(i)], 'r--','LineWidth',1);
%         plot(PlotScore(i), y(i), 'ro','MarkerFaceColor','r','MarkerSize',6);
%     elseif Plotp(i) <= 0.05
%         plot([0 PlotScore(i)], [y(i) y(i)], 'k--','LineWidth',1);
%         plot(PlotScore(i), y(i), 'ko','MarkerFaceColor','k','MarkerSize',6);
%     else
%         plot([0 PlotScore(i)], [y(i) y(i)], 'k--','LineWidth',1);
%         plot(PlotScore(i), y(i), 'ko','MarkerFaceColor','k','MarkerSize',6);
%     end
% end
% plot([0 0],[0 numel(PlotScore)+1],'-','color',[0.7 0.7 0.7],'LineWidth',1);
% set(gca,'xLim',[-0.8 0.8],'xTick',[-0.6 0 0.6],'FontName', 'Arial');
% set(gca,'yLim',[0.5 numel(PlotScore)+0.5],'yTick',[1:numel(PlotScore)],'yTickLabel',cellLabel(y), 'FontName', 'Arial');
% set(gca, 'Box', 'off');
% ax = gca;
% ax.YColor = 'k'; 
% ax.XColor = 'k'; % 保持刻度数字
% ax.TickDir = 'out';





%% plotting scatter
CorrType = 'Pearson';  % Pearson Spearman
figure('Color',[1 1 1],'Position',[0 0 400 600],'Units','pixels');
% which_maps = [1 18 25];
which_maps = [2 18 19];

% CBF-score
Data_name = ['CBF_PC1_Val_L' num2str(num_layer)];
origx = load([Group_dir Data_name UniOrBil '.1D']);
origData1_null = load([Group_dir Data_name UniOrBil nullType '.1D']);
GeneMask = GeneMask_SNR0;
Sample_Mask = Sample_Mask_SNR0;

r = [];
spinp = [];
for i = 1:numel(which_maps)
    y = Cell_dis_m2(:,which_maps(i));
    x = origx(Sample_Mask);
    y = y(GeneMask);
    Data1_null = origData1_null(Sample_Mask,:);

    c = corr(y,x,'Type',CorrType);
    cn = corr(Data1_null,y,'Type',CorrType);
    p = (1+sum(abs(c)<abs(cn)))/20001;
    r = [r; c];
    spinp = [spinp; p];

    subplot(numel(which_maps),2,i*2-1); hold on;
    subfun_plot_scatter_regression( x , y , [-10 8] , [floor(min(y*10))/10 ceil(max(y*10))/10]);
end
r_CBF = r;
p_CBF = spinp;


% CCSI
Data_name = ['CBF_Mean_PScs_L' num2str(num_layer)];
origx = load([Group_dir Data_name UniOrBil '.1D']);
origData1_null = load([Group_dir Data_name UniOrBil nullType '.1D']);
GeneMask = GeneMask_SNR4;
Sample_Mask = Sample_Mask_SNR4;

r = [];
spinp = [];
for i = 1:numel(which_maps)
    y = Cell_dis_m2(:,which_maps(i));
    x = origx(Sample_Mask);
    y = y(GeneMask);
    Data1_null = origData1_null(Sample_Mask,:);

    c = corr(y,x,'Type',CorrType);
    cn = corr(Data1_null,y,'Type',CorrType);
    p = (1+sum(abs(c)<abs(cn)))/20001;
    r = [r; c];
    spinp = [spinp; p];

    subplot(numel(which_maps),2,i*2); hold on;
    subfun_plot_scatter_regression( x , y , [-0.2 1.2] , [floor(min(y*10))/10 ceil(max(y*10))/10]);
end
r_CCSI = r;
p_CCSI = spinp;


%% write to map plot by R
OutNames = {'Endo_uni.csv', 'Per_uni.csv',...
            'In6a_uni.csv', 'In6b_uni.csv',...
            'Oligo_uni.csv'};
OutDatas = {zscore(Cell_dis_m2(:,1)), zscore(Cell_dis_m2(:,2)),...
            zscore(Cell_dis_m2(:,18)), zscore(Cell_dis_m2(:,19)),...
            zscore(Cell_dis_m2(:,25))};
Datalabel = Sample_NumLable;

% for i = 1:numel(OutNames)
%     tmp = subfun_write_to_MapPlot_byR(OutDatas{i},Datalabel,atlas_info,[Group_dir OutNames{i}]);
% end



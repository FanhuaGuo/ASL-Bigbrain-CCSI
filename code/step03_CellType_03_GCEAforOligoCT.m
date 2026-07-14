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
CT_categories3_orig = CT_categories3;


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


%% get Oligo-cate without any overlap
ifDeleteOverlap = 0;
if ifDeleteOverlap
    CT_categories3 = CT_categories3_orig(3:8,:);
    tmpGene = zeros(numel(GeneID),1);
    for i = 1:size(CT_categories3,1)
        tmpGene(CT_categories3.gene_label{i}) = tmpGene(CT_categories3.gene_label{i})+1;
    end
    for i = 1:size(CT_categories3,1)
        for j = 1:numel(CT_categories3.gene_label{i})
            if tmpGene(CT_categories3.gene_label{i}(j)) > 1
                CT_categories3.gene_label{i}(j) = 0;
            end
        end
        CT_categories3.gene_label{i}(CT_categories3.gene_label{i}==0) = [];
    end
else
    CT_categories3 = CT_categories3_orig;
end


%% ======================= Method1 median correlation's r-fisherZ ========================
%  =======================================================================================
%% common parameters
UniOrBil = '_uni';
BH_alpha = 0.05;
Strategy = 2; % 1:direct mean;  2:direct median
cate = CT_categories3.gene_label;
cateLabel = CT_categories3.Cell;
CorrType = 'Pearson';  % Pearson Spearman
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


%% update FDR
p = [GCEA_CBFscore_1.spin_p; GCEA_CCSI_1.spin_p; GCEA_CCSIcon_1.spin_p];
FDR = fdr_BH(p,0.5);
GCEA_CBFscore_1.FDRtog = FDR(1:end/3)';
GCEA_CCSI_1.FDRtog = FDR(end/3+1:2*end/3)';
GCEA_CCSIcon_1.FDRtog = FDR(2*end/3+1:end)';



%% plotting v3
% manual parameters
DataPlotcbf = GCEA_CBFscore_1;
DataPlotccsi = GCEA_CCSI_1;
DataPlotccsicon = GCEA_CCSIcon_1;
cellLabel = GCEA_CCSI_1.Cell;


% auto parameters
PlotScore = [DataPlotcbf.GCEAscore DataPlotccsi.GCEAscore DataPlotccsicon.GCEAscore];
PlotFDR = [DataPlotcbf.FDRtog DataPlotccsi.FDRtog DataPlotccsicon.FDRtog];
if ifDeleteOverlap
    PlotScore = PlotScore(end:-1:1,:);
    PlotFDR = PlotFDR(end:-1:1,:);
    cellLabel = cellLabel(end:-1:1);
else
    PlotScore = PlotScore([8:-1:3 9 2 1],:);
    PlotFDR = PlotFDR([8:-1:3 9 2 1],:);
    cellLabel = cellLabel([8:-1:3 9 2 1]);
end
tCellNum = size(PlotScore,1);


% plotting
figure('Color',[1 1 1],'Position',[0 0 400 400],'Units','pixels');
hold on;
shift = 0.5;
for i = 1:tCellNum
    if PlotFDR(i,1) <= 0.01
        plot([0 PlotScore(i,1)], [i*2-1+shift i*2-1+shift], 'r--','LineWidth',0.25);
        plot(PlotScore(i,1), i*2-1+shift, 'ro','MarkerFaceColor','r','MarkerSize',6);
    else
        plot([0 PlotScore(i,1)], [i*2-1+shift i*2-1+shift], 'k--','LineWidth',0.25);
        plot(PlotScore(i,1), i*2-1+shift, 'ko','MarkerFaceColor','k','MarkerSize',6);
    end

    if PlotFDR(i,2) <= 0.01
        plot([0 PlotScore(i,2)], [i*2-1 i*2-1], 'r--','LineWidth',0.25);
        plot(PlotScore(i,2), i*2-1, 'r^','MarkerFaceColor','r','MarkerSize',6);
    else
        plot([0 PlotScore(i,2)], [i*2-1 i*2-1], 'k--','LineWidth',0.25);
        plot(PlotScore(i,2), i*2-1, 'k^','MarkerFaceColor','k','MarkerSize',6);
    end
    
    if PlotFDR(i,3) <= 0.01
        plot([0 PlotScore(i,3)], [i*2-1-shift i*2-1-shift], 'r--','LineWidth',0.25);
        plot(PlotScore(i,3), i*2-1-shift, 'rv','MarkerFaceColor','r','MarkerSize',6);
    else
        plot([0 PlotScore(i,3)], [i*2-1-shift i*2-1-shift], 'k--','LineWidth',0.25);
        plot(PlotScore(i,3), i*2-1-shift, 'kv','MarkerFaceColor','k','MarkerSize',6);
    end
end
plot([0 0],[0 2*tCellNum+1],'-','color',[0.7 0.7 0.7],'LineWidth',0.25);
set(gca,'xLim',[-0.4 0.4],'xTick',[-0.3 0 0.3],'FontName', 'Arial');
set(gca,'yLim',[0 2*tCellNum],'yTick',[1:2:2*tCellNum],'yTickLabel',cellLabel, 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';


%% plotting v2 (old version)
% % manual parameters
% DataPlotcbf = GCEA_CBFscore_1;
% DataPlotccsi = GCEA_CCSI_1;
% cellLabel = GCEA_CCSI_1.Cell;
% 
% 
% % auto parameters
% PlotScore = [DataPlotcbf.GCEAscore DataPlotccsi.GCEAscore];
% Plotp = [DataPlotcbf.spin_p DataPlotccsi.spin_p];
% if ifDeleteOverlap
%     PlotScore = PlotScore(end:-1:1,:);
%     Plotp = Plotp(end:-1:1,:);
%     cellLabel = cellLabel(end:-1:1);
% else
%     PlotScore = PlotScore([8:-1:3 9 2 1],:);
%     Plotp = Plotp([8:-1:3 9 2 1],:);
%     cellLabel = cellLabel([8:-1:3 9 2 1]);
% end
% PlotFDR = [fdr_BH(Plotp(:,1),0.05)' fdr_BH(Plotp(:,2),0.05)'];
% tCellNum = size(PlotScore,1);
% 
% 
% % plotting
% figure('Color',[1 1 1],'Position',[0 0 400 400],'Units','pixels');
% hold on;
% for i = 1:tCellNum
%     if PlotFDR(i,1) <= 0.05
%         plot([0 PlotScore(i,1)], [i+0.15 i+0.15], 'r--','LineWidth',1);
%         plot(PlotScore(i,1), i+0.15, 'ro','MarkerFaceColor','r','MarkerSize',6);
%     elseif Plotp(i,1) <= 0.05
%         plot([0 PlotScore(i,1)], [i+0.15 i+0.15], 'k--','LineWidth',1);
%         plot(PlotScore(i,1), i+0.15, 'ko','MarkerFaceColor','k','MarkerSize',6);
%     else
%         plot([0 PlotScore(i,1)], [i+0.15 i+0.15], 'k--','LineWidth',1);
%         plot(PlotScore(i,1), i+0.15, 'ko','MarkerFaceColor','k','MarkerSize',6);
%     end
% 
%     if PlotFDR(i,2) <= 0.05
%         plot([0 PlotScore(i,2)], [i-0.15 i-0.15], 'r--','LineWidth',1);
%         plot(PlotScore(i,2), i-0.15, 'r^','MarkerFaceColor','r','MarkerSize',6);
%     elseif Plotp(i,2) <= 0.05
%         plot([0 PlotScore(i,2)], [i-0.15 i-0.15], 'k--','LineWidth',1);
%         plot(PlotScore(i,2), i-0.15, 'k^','MarkerFaceColor','k','MarkerSize',6);
%     else
%         plot([0 PlotScore(i,2)], [i-0.15 i-0.15], 'k--','LineWidth',1);
%         plot(PlotScore(i,2), i-0.15, 'k^','MarkerFaceColor','k','MarkerSize',6);
%     end
% end
% plot([0 0],[0 tCellNum+1],'-','color',[0.7 0.7 0.7],'LineWidth',1);
% set(gca,'xLim',[-0.4 0.4],'xTick',[-0.3 0 0.3],'FontName', 'Arial');
% set(gca,'yLim',[0.5 tCellNum+0.5],'yTick',[1:tCellNum],'yTickLabel',cellLabel, 'FontName', 'Arial');
% set(gca, 'Box', 'off');
% ax = gca;
% ax.YColor = 'k'; 
% ax.XColor = 'k'; % 保持刻度数字
% ax.TickDir = 'out';




%% ======================= Method2 median Expression maps ================================
%  =======================================================================================
%% compute cell distribution by median GeneExp
% common parameters
cate = CT_categories3.gene_label;
cateLabel = CT_categories3.Cell;
UniOrBil = '_uni';
BH_alpha = 0.05;
CorrType = 'Pearson';  % Pearson Spearman

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



%% update FDR
p = [GCEA_CBFscore_2.spin_p; GCEA_CCSI_2.spin_p; GCEA_CCSIcon_2.spin_p];
FDR = fdr_BH(p,0.5);
GCEA_CBFscore_2.FDRtog = FDR(1:end/3)';
GCEA_CCSI_2.FDRtog = FDR(end/3+1:2*end/3)';
GCEA_CCSIcon_2.FDRtog = FDR(2*end/3+1:end)';



%% plotting v3
% manual parameters
DataPlotcbf = GCEA_CBFscore_2;
DataPlotccsi = GCEA_CCSI_2;
DataPlotccsicon = GCEA_CCSIcon_2;
cellLabel = GCEA_CCSI_2.Cell;


% auto parameters
PlotScore = [DataPlotcbf.GCEAscore DataPlotccsi.GCEAscore DataPlotccsicon.GCEAscore];
PlotFDR = [DataPlotcbf.FDRtog DataPlotccsi.FDRtog DataPlotccsicon.FDRtog];
if ifDeleteOverlap
    PlotScore = PlotScore(end:-1:1,:);
    PlotFDR = PlotFDR(end:-1:1,:);
    cellLabel = cellLabel(end:-1:1);
else
    PlotScore = PlotScore([8:-1:3 9 2 1],:);
    PlotFDR = PlotFDR([8:-1:3 9 2 1],:);
    cellLabel = cellLabel([8:-1:3 9 2 1]);
end
tCellNum = size(PlotScore,1);


% plotting
figure('Color',[1 1 1],'Position',[0 0 400 400],'Units','pixels');
hold on;
shift = 0.5;
for i = 1:tCellNum
    if PlotFDR(i,1) <= 0.05
        plot([0 PlotScore(i,1)], [i*2-1+shift i*2-1+shift], 'r--','LineWidth',0.25);
        plot(PlotScore(i,1), i*2-1+shift, 'ro','MarkerFaceColor','r','MarkerSize',6);
    else
        plot([0 PlotScore(i,1)], [i*2-1+shift i*2-1+shift], 'k--','LineWidth',0.25);
        plot(PlotScore(i,1), i*2-1+shift, 'ko','MarkerFaceColor','k','MarkerSize',6);
    end

    if PlotFDR(i,2) <= 0.05
        plot([0 PlotScore(i,2)], [i*2-1 i*2-1], 'r--','LineWidth',0.25);
        plot(PlotScore(i,2), i*2-1, 'r^','MarkerFaceColor','r','MarkerSize',6);
    else
        plot([0 PlotScore(i,2)], [i*2-1 i*2-1], 'k--','LineWidth',0.25);
        plot(PlotScore(i,2), i*2-1, 'k^','MarkerFaceColor','k','MarkerSize',6);
    end
    
    if PlotFDR(i,3) <= 0.05
        plot([0 PlotScore(i,3)], [i*2-1-shift i*2-1-shift], 'r--','LineWidth',0.25);
        plot(PlotScore(i,3), i*2-1-shift, 'rv','MarkerFaceColor','r','MarkerSize',6);
    else
        plot([0 PlotScore(i,3)], [i*2-1-shift i*2-1-shift], 'k--','LineWidth',0.25);
        plot(PlotScore(i,3), i*2-1-shift, 'kv','MarkerFaceColor','k','MarkerSize',6);
    end
end
plot([0 0],[0 2*tCellNum+1],'-','color',[0.7 0.7 0.7],'LineWidth',0.25);
set(gca,'xLim',[-0.6 0.6],'xTick',[-0.4 0 0.4],'FontName', 'Arial');
set(gca,'yLim',[0 2*tCellNum],'yTick',[1:2:2*tCellNum],'yTickLabel',cellLabel, 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';




%% plotting v2 (old version)
% % manual parameters
% DataPlotcbf = GCEA_CBFscore_2;
% DataPlotccsi = GCEA_CCSI_2;
% cellLabel = GCEA_CCSI_2.Cell;
% 
% 
% % auto parameters
% PlotScore = [DataPlotcbf.GCEAscore DataPlotccsi.GCEAscore];
% Plotp = [DataPlotcbf.spin_p DataPlotccsi.spin_p];
% PlotScore = PlotScore([8:-1:3 9 2 1],:);
% Plotp = Plotp([8:-1:3 9 2 1],:);
% PlotFDR = [fdr_BH(Plotp(:,1),0.05)' fdr_BH(Plotp(:,2),0.05)'];
% tCellNum = size(PlotScore,1);
% 
% 
% % plotting
% figure('Color',[1 1 1],'Position',[0 0 400 400],'Units','pixels');
% hold on;
% for i = 1:tCellNum
%     if PlotFDR(i,1) <= 0.05
%         plot([0 PlotScore(i,1)], [i+0.15 i+0.15], 'r--','LineWidth',1);
%         plot(PlotScore(i,1), i+0.15, 'ro','MarkerFaceColor','r','MarkerSize',6);
%     elseif Plotp(i,1) <= 0.05
%         plot([0 PlotScore(i,1)], [i+0.15 i+0.15], 'k--','LineWidth',1);
%         plot(PlotScore(i,1), i+0.15, 'ko','MarkerFaceColor','k','MarkerSize',6);
%     else
%         plot([0 PlotScore(i,1)], [i+0.15 i+0.15], 'k--','LineWidth',1);
%         plot(PlotScore(i,1), i+0.15, 'ko','MarkerFaceColor','k','MarkerSize',6);
%     end
% 
%     if PlotFDR(i,2) <= 0.05
%         plot([0 PlotScore(i,2)], [i-0.15 i-0.15], 'r--','LineWidth',1);
%         plot(PlotScore(i,2), i-0.15, 'r^','MarkerFaceColor','r','MarkerSize',6);
%     elseif Plotp(i,2) <= 0.05
%         plot([0 PlotScore(i,2)], [i-0.15 i-0.15], 'k--','LineWidth',1);
%         plot(PlotScore(i,2), i-0.15, 'k^','MarkerFaceColor','k','MarkerSize',6);
%     else
%         plot([0 PlotScore(i,2)], [i-0.15 i-0.15], 'k--','LineWidth',1);
%         plot(PlotScore(i,2), i-0.15, 'k^','MarkerFaceColor','k','MarkerSize',6);
%     end
% end
% plot([0 0],[0 tCellNum+1],'-','color',[0.7 0.7 0.7],'LineWidth',1);
% set(gca,'xLim',[-0.6 0.6],'xTick',[-0.4 0 0.4],'FontName', 'Arial');
% set(gca,'yLim',[0.5 tCellNum+0.5],'yTick',[1:tCellNum],'yTickLabel',cellLabel([8:-1:3 9 2 1]), 'FontName', 'Arial');
% set(gca, 'Box', 'off');
% ax = gca;
% ax.YColor = 'k'; 
% ax.XColor = 'k'; % 保持刻度数字
% ax.TickDir = 'out';




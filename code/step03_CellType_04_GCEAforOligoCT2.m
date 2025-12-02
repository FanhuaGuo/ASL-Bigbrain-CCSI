clc
clear all
close all
set(0,'defaultfigurecolor',[1 1 1]);% set figure background color
addpath(genpath('/Users/guofanhua/Desktop/gfh/tools/MatlabTools/BrainSpace-latest/matlab'))

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
ifDeleteOverlap = 1;
if ifDeleteOverlap
    CT_categories3 = CT_categories3_orig(3:8,:);
    tmpGene = zeros(numel(GeneID),6);
    clc
    
    % record gene id
    fprintf('Gene set number(Oligo 1->6):\n');
    for i = 1:size(CT_categories3,1)
        tmpGene(CT_categories3.gene_label{i},i) = 1;
        fprintf([num2str(numel(CT_categories3.gene_label{i})) ' ']);
    end
    
    % check overlap number
    Gene_Overlap = zeros(6);
    for i = 1:6
        for j = 1:6
            Gene_Overlap(i,j) = sum(tmpGene(:,i)==1 & tmpGene(:,j)==1);
        end
    end
    
    % delete overlap
    fprintf('\n');
    fprintf('Gene set number(Oligo 1->6) without overlap:\n');
    for i = 1:size(CT_categories3,1)
        for j = 1:numel(CT_categories3.gene_label{i})
            if sum(tmpGene(CT_categories3.gene_label{i}(j),:)) > 1
                CT_categories3.gene_label{i}(j) = 0;
            end
        end
        CT_categories3.gene_label{i}(CT_categories3.gene_label{i}==0) = [];
        fprintf([num2str(numel(CT_categories3.gene_label{i})) ' ']);
    end
    fprintf('\n');
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


%% compute CBF CellType1-GSEA by median correlation's r-fisherZ
DataName = 'CBF_PC1_Val';
GeneMask = GeneMask_SNR0;
Sample_Mask = Sample_Mask_SNR0;


% read data
RealData_path = [Group_dir DataName UniOrBil '.1D'];
NullData_path = [Group_dir DataName UniOrBil '_null.1D'];
GeneData = GeneExpression_Orig(:,GeneMask);
MapData = load(RealData_path);
MapData_null = load(NullData_path);
RepeatNum = size(MapData_null,2);


% compute GSEA
r = corr(MapData(Sample_Mask),GeneData','Type',CorrType);
r_Fz = subfun_fisher_z(r);
rn = corr(MapData_null(Sample_Mask,:),GeneData','Type',CorrType);
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
DataName = 'CBF_Mean_PScs';
GeneMask = GeneMask_SNR4;
Sample_Mask = Sample_Mask_SNR4;


% read data
RealData_path = [Group_dir DataName UniOrBil '.1D'];
NullData_path = [Group_dir DataName UniOrBil '_null.1D'];
GeneData = GeneExpression_Orig(:,GeneMask);
MapData = load(RealData_path);
MapData_null = load(NullData_path);
RepeatNum = size(MapData_null,2);


% compute GSEA
r = corr(MapData(Sample_Mask),GeneData','Type',CorrType);
r_Fz = subfun_fisher_z(r);
rn = corr(MapData_null(Sample_Mask,:),GeneData','Type',CorrType);
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

      

%% plotting cbf
% manual parameters
DataPlot = GCEA_CBFscore_1;
cellLabel = GCEA_CBFscore_1.Cell;


% auto parameters
PlotScore = DataPlot.GCEAscore;
Plotp = DataPlot.spin_p;
PlotScore = PlotScore([3 4 6 1 5 2]);
Plotp = Plotp([3 4 6 1 5 2]);
cellLabel = cellLabel([3 4 6 1 5 2]);

PlotFDR = fdr_BH(Plotp,0.05)';
tCellNum = size(PlotScore,1);


% plotting
figure('Color',[1 1 1],'Position',[0 0 300 200],'Units','pixels');
hold on;
b = bar(1:2, PlotScore(1:2), 'FaceColor', 'white', 'EdgeColor', [255, 31, 201]/255, 'LineWidth', 1);
b = bar(3, PlotScore(3), 'FaceColor', 'white', 'EdgeColor', [0, 191, 86]/255, 'LineWidth', 1);
b = bar(4:5, PlotScore(4:5), 'FaceColor', 'white', 'EdgeColor', [19, 0, 170]/255, 'LineWidth', 1);
b = bar(6, PlotScore(6), 'FaceColor', 'white', 'EdgeColor', [0, 197, 212]/255, 'LineWidth', 1);
set(gca,'xLim',[0.5 tCellNum+0.5],'xTick',[1:tCellNum],'xTickLabel',cellLabel, 'FontName', 'Arial');
set(gca, 'ylim', [-0.2 0.6], 'yTick', -0.2:0.2:0.6, 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';



%% plotting ccsi
% manual parameters
DataPlot = GCEA_CCSI_1;
cellLabel = GCEA_CCSI_1.Cell;


% auto parameters
PlotScore = DataPlot.GCEAscore;
Plotp = DataPlot.spin_p;
PlotScore = PlotScore([3 4 6 1 5 2]);
Plotp = Plotp([3 4 6 1 5 2]);
cellLabel = cellLabel([3 4 6 1 5 2]);

PlotFDR = fdr_BH(Plotp,0.05)';
tCellNum = size(PlotScore,1);


% plotting
figure('Color',[1 1 1],'Position',[0 0 300 200],'Units','pixels');
hold on;
b = bar(1:2, PlotScore(1:2), 'FaceColor', 'white', 'EdgeColor', [255, 31, 201]/255, 'LineWidth', 1);
b = bar(3, PlotScore(3), 'FaceColor', 'white', 'EdgeColor', [0, 191, 86]/255, 'LineWidth', 1);
b = bar(4:5, PlotScore(4:5), 'FaceColor', 'white', 'EdgeColor', [19, 0, 170]/255, 'LineWidth', 1);
b = bar(6, PlotScore(6), 'FaceColor', 'white', 'EdgeColor', [0, 197, 212]/255, 'LineWidth', 1);
set(gca,'xLim',[0.5 tCellNum+0.5],'xTick',[1:tCellNum],'xTickLabel',cellLabel, 'FontName', 'Arial');
set(gca, 'ylim', [-0.2 0.4], 'yTick', -0.2:0.2:0.6, 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';




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
DataName = 'CBF_PC1_Val';
GeneMask = GeneMask_SNR0;
Sample_Mask = Sample_Mask_SNR0;


% read data
RealData_path = [Group_dir DataName UniOrBil '.1D'];
NullData_path = [Group_dir DataName UniOrBil '_null.1D'];
MapData = load(RealData_path);
MapData_null = load(NullData_path);
RepeatNum = size(MapData_null,2);


% compute GSEA
r = corr(MapData(Sample_Mask),Cell_dis_m2(GeneMask,:),'Type',CorrType);
r_Fz = subfun_fisher_z(r);
rn = corr(MapData_null(Sample_Mask,:),Cell_dis_m2(GeneMask,:),'Type',CorrType);
rn_Fz = subfun_fisher_z(rn);
spin_p = (1+sum(abs(r_Fz) < abs(rn_Fz)))/(1+RepeatNum);
BHp = fdr_BH(spin_p,BH_alpha)';


% disposal result
GCEAresults = table(cateLabel, r_Fz', spin_p', BHp,...
                    'VariableNames',{'Cell', 'GCEAscore', 'spin_p', 'FDR'});

GCEA_CBFscore_2 = GCEAresults;
                
      


%% compute CCSI CellType1-GSEA by median GeneExp as Cell-distribution
DataName = 'CBF_Mean_PScs';
GeneMask = GeneMask_SNR4;
Sample_Mask = Sample_Mask_SNR4;


% read data
RealData_path = [Group_dir DataName UniOrBil '.1D'];
NullData_path = [Group_dir DataName UniOrBil '_null.1D'];
MapData = load(RealData_path);
MapData_null = load(NullData_path);
RepeatNum = size(MapData_null,2);


% compute GSEA
r = corr(MapData(Sample_Mask),Cell_dis_m2(GeneMask,:),'Type',CorrType);
r_Fz = subfun_fisher_z(r);
rn = corr(MapData_null(Sample_Mask,:),Cell_dis_m2(GeneMask,:),'Type',CorrType);
rn_Fz = subfun_fisher_z(rn);
spin_p = (1+sum(abs(r_Fz) < abs(rn_Fz)))/(1+RepeatNum);
BHp = fdr_BH(spin_p,BH_alpha)';


% disposal result
GCEAresults = table(cateLabel, r_Fz', spin_p', BHp,...
                    'VariableNames',{'Cell', 'GCEAscore', 'spin_p', 'FDR'});

GCEA_CCSI_2 = GCEAresults;



%% plotting cbf
% manual parameters
DataPlot = GCEA_CBFscore_2;
cellLabel = GCEA_CBFscore_2.Cell;


% auto parameters
PlotScore = DataPlot.GCEAscore;
Plotp = DataPlot.spin_p;
PlotScore = PlotScore([3 4 6 1 5 2]);
Plotp = Plotp([3 4 6 1 5 2]);
cellLabel = cellLabel([3 4 6 1 5 2]);

PlotFDR = fdr_BH(Plotp,0.05)';
tCellNum = size(PlotScore,1);


% plotting
figure('Color',[1 1 1],'Position',[0 0 300 200],'Units','pixels');
hold on;
b = bar(1:2, PlotScore(1:2), 'FaceColor', 'white', 'EdgeColor', [0, 200, 180]/255, 'LineWidth', 1);
b = bar(3, PlotScore(3), 'FaceColor', 'white', 'EdgeColor', [242, 142, 43]/255, 'LineWidth', 1);
b = bar(4:5, PlotScore(4:5), 'FaceColor', 'white', 'EdgeColor', [89, 161, 79]/255, 'LineWidth', 1);
b = bar(6, PlotScore(6), 'FaceColor', 'white', 'EdgeColor', [225, 87, 89]/255, 'LineWidth', 1);
set(gca,'xLim',[0.5 tCellNum+0.5],'xTick',[1:tCellNum],'xTickLabel',cellLabel, 'FontName', 'Arial');
set(gca, 'ylim', [-0.2 0.6], 'yTick', -0.2:0.2:0.6, 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';



%% plotting ccsi
% manual parameters
DataPlot = GCEA_CCSI_2;
cellLabel = GCEA_CCSI_2.Cell;


% auto parameters
PlotScore = DataPlot.GCEAscore;
Plotp = DataPlot.spin_p;
PlotScore = PlotScore([3 4 6 1 5 2]);
Plotp = Plotp([3 4 6 1 5 2]);
cellLabel = cellLabel([3 4 6 1 5 2]);

PlotFDR = fdr_BH(Plotp,0.05)';
tCellNum = size(PlotScore,1);


% plotting
figure('Color',[1 1 1],'Position',[0 0 300 200],'Units','pixels');
hold on;
b = bar(1:2, PlotScore(1:2), 'FaceColor', 'white', 'EdgeColor', [0, 200, 180]/255, 'LineWidth', 1);
b = bar(3, PlotScore(3), 'FaceColor', 'white', 'EdgeColor', [242, 142, 43]/255, 'LineWidth', 1);
b = bar(4:5, PlotScore(4:5), 'FaceColor', 'white', 'EdgeColor', [89, 161, 79]/255, 'LineWidth', 1);
b = bar(6, PlotScore(6), 'FaceColor', 'white', 'EdgeColor', [225, 87, 89]/255, 'LineWidth', 1);
set(gca,'xLim',[0.5 tCellNum+0.5],'xTick',[1:tCellNum],'xTickLabel',cellLabel, 'FontName', 'Arial');
set(gca, 'ylim', [-0.2 0.6], 'yTick', -0.2:0.2:0.6, 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';





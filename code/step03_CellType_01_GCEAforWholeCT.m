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


%% compute CBF CellType1-GSEA by median correlation's r-fisherZ
DataName = 'CBF_PC3_Val';
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



%% plotting CCSI
% manual parameters
DataPlot = GCEA_CCSI_2;


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





%% plotting scatter
CorrType = 'Pearson';  % Pearson Spearman
figure('Color',[1 1 1],'Position',[0 0 400 600],'Units','pixels');
% which_maps = [1 18 25];
which_maps = [2 18 19];

% CBF-score
Data_name = 'CBF_PC1_Val_uni';  % 'CBF_PC1_Val_bil', 'CBF_Mean_PScs_bil'
origx = load([Group_dir Data_name '.1D']);
origData1_null = load([Group_dir Data_name '_null.1D']);
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
Data_name = 'CBF_Mean_PScs_uni';  % 'CBF_PC1_Val_bil', 'CBF_Mean_PScs_bil'
origx = load([Group_dir Data_name '.1D']);
origData1_null = load([Group_dir Data_name '_null.1D']);
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
    subfun_plot_scatter_regression( x , y , [0 1.2] , [floor(min(y*10))/10 ceil(max(y*10))/10]);
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



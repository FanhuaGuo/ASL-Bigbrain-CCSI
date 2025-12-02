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
GCEAresults = table(GO_categories.go_id, GO_categories.go_func, cate_rZ', spin_p', BHp, GeneNum_GO,...
                    'VariableNames',{'go_id', 'go_func', 'GCEAscore', 'spin_p', 'FDR', 'GeneNum'});

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
GCEAresults = table(GO_categories.go_id, GO_categories.go_func, cate_rZ', spin_p', BHp, GeneNum_GO,...
                    'VariableNames',{'go_id', 'go_func', 'GCEAscore', 'spin_p', 'FDR', 'GeneNum'});

GCEA_CCSI_1 = GCEAresults;



%% plotting CBF
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




%% ======================= Method2 median Expression maps ================================
%  =======================================================================================
% %% compute cell distribution by median GeneExp
% % common parameters
% UniOrBil = '_uni';
% BH_alpha = 0.05;
% cate = GO_categories.gene_label;
% 
% 
% 
% % compute
% GO_dis_m2 = [];
% GeneData = GeneExpression_Orig;
% for i = 1:numel(cate)
%     GO_dis_m2 = [GO_dis_m2; median(GeneExpression_Orig(cate{i},:),1)];
% end
% GO_dis_m2 = GO_dis_m2';
% 
% 
% %% compute CBF CellType1-GSEA by median GeneExp as Cell-distribution
% DataName = 'CBF_PC1_Val';
% GeneMask = GeneMask_SNR0;
% Sample_Mask = Sample_Mask_SNR0;
% 
% 
% % read data
% RealData_path = [Group_dir DataName UniOrBil '.1D'];
% NullData_path = [Group_dir DataName UniOrBil '_null.1D'];
% MapData = load(RealData_path);
% MapData_null = load(NullData_path);
% RepeatNum = size(MapData_null,2);
% 
% 
% % compute GSEA
% r = corr(MapData(Sample_Mask),GO_dis_m2(GeneMask,:),'Type',CorrType);
% r_Fz = subfun_fisher_z(r);
% rn = corr(MapData_null(Sample_Mask,:),GO_dis_m2(GeneMask,:),'Type',CorrType);
% rn_Fz = subfun_fisher_z(rn);
% spin_p = (1+sum(abs(r_Fz) < abs(rn_Fz)))/(1+RepeatNum);
% BHp = fdr_BH(spin_p,BH_alpha)';
% 
% 
% % disposal result
% GCEAresults = table(GO_categories.go_id, GO_categories.go_func, cate_rZ', spin_p', BHp, GeneNum_GO,...
%                     'VariableNames',{'go_id', 'go_func', 'GCEAscore', 'spin_p', 'FDR', 'GeneNum'});
% 
% GCEA_CBFscore_2 = GCEAresults;
%                 
%       
% 
% 
% %% compute CCSI CellType1-GSEA by median GeneExp as Cell-distribution
% DataName = 'CBF_Mean_PScs';
% GeneMask = GeneMask_SNR4;
% Sample_Mask = Sample_Mask_SNR4;
% 
% 
% % read data
% RealData_path = [Group_dir DataName UniOrBil '.1D'];
% NullData_path = [Group_dir DataName UniOrBil '_null.1D'];
% MapData = load(RealData_path);
% MapData_null = load(NullData_path);
% RepeatNum = size(MapData_null,2);
% 
% 
% % compute GSEA
% r = corr(MapData(Sample_Mask),GO_dis_m2(GeneMask,:),'Type',CorrType);
% r_Fz = subfun_fisher_z(r);
% rn = corr(MapData_null(Sample_Mask,:),GO_dis_m2(GeneMask,:),'Type',CorrType);
% rn_Fz = subfun_fisher_z(rn);
% spin_p = (1+sum(abs(r_Fz) < abs(rn_Fz)))/(1+RepeatNum);
% BHp = fdr_BH(spin_p,BH_alpha)';
% 
% 
% % disposal result
% GCEAresults = table(GO_categories.go_id, GO_categories.go_func, cate_rZ', spin_p', BHp, GeneNum_GO,...
%                     'VariableNames',{'go_id', 'go_func', 'GCEAscore', 'spin_p', 'FDR', 'GeneNum'});
% 
% GCEA_CCSI_2 = GCEAresults;
% 
% 
% 
% %% plotting CBF
% Plotdata = GCEA_CBFscore_2;
% 
% Plotdata = Plotdata(Plotdata.FDR <= 0.05,:);
% Plotdata = sortrows(Plotdata,'spin_p','ascend');
% Plotscore = Plotdata.GCEAscore;
% PlotSize = Plotdata.GeneNum;
% PlotSize(PlotSize>100) = 100;
% PlotFDR = [];
% for i = 1:size(Plotdata,1)
%     if Plotdata.spin_p(i) == 1/20001
%         PlotFDR = [PlotFDR; 3.2];
%     else
%         PlotFDR = [PlotFDR; -log10(Plotdata.FDR(i))];
%     end
% end
% 
% 
% % plotting
% figure('Color',[1 1 1],'Position',[0 0 400 250],'Units','pixels');
% hold on;
% scatter(PlotFDR, Plotscore, PlotSize*5, 'fill',...
%         'MarkerFaceColor', [0.7 0.7 0.7]);
% scatter([3.2; 3.2; 3.2], [-0.3; -0.2; -0.1], [20; 60; 100]*5, 'fill',...
%         'MarkerFaceColor', [1 0 0]);
% set(gca,'xLim',[1 3.5],'xTick',[-log10(0.05) -log10(0.01) -log10(0.001) 3.2], 'xTickLabel', {'0.05','0.01','0.001','spin 0'},'FontName', 'Arial');
% set(gca,'yLim',[-0.4 0.4],'yTick',[-0.4 -0.3 0 0.3 0.4], 'FontName', 'Arial');
% set(gca, 'Box', 'off');
% ax = gca;
% ax.YColor = 'k'; 
% ax.XColor = 'k'; % 保持刻度数字
% ax.TickDir = 'out';
% 
% 
% 
% 
% 
% %% plotting CCSI
% Plotdata = GCEA_CCSI_2;
% 
% Plotdata = Plotdata(Plotdata.FDR <= 0.05,:);
% Plotdata = sortrows(Plotdata,'spin_p','ascend');
% Plotscore = Plotdata.GCEAscore;
% PlotSize = Plotdata.GeneNum;
% PlotSize(PlotSize>100) = 100;
% PlotFDR = [];
% for i = 1:size(Plotdata,1)
%     if Plotdata.spin_p(i) == 1/20001
%         PlotFDR = [PlotFDR; 3.2];
%     else
%         PlotFDR = [PlotFDR; -log10(Plotdata.FDR(i))];
%     end
% end
% 
% 
% % plotting
% figure('Color',[1 1 1],'Position',[0 0 400 250],'Units','pixels');
% hold on;
% scatter(PlotFDR, Plotscore, PlotSize*5, 'fill',...
%         'MarkerFaceColor', [0.7 0.7 0.7]);
% scatter([3.2; 3.2; 3.2], [-0.3; -0.2; -0.1], [20; 60; 100]*5, 'fill',...
%         'MarkerFaceColor', [1 0 0]);
% set(gca,'xLim',[1 3.5],'xTick',[-log10(0.05) -log10(0.01) -log10(0.001) 3.2], 'xTickLabel', {'0.05','0.01','0.001','spin 0'},'FontName', 'Arial');
% set(gca,'yLim',[-0.4 0.4],'yTick',[-0.4 -0.3 0 0.3 0.4], 'FontName', 'Arial');
% set(gca, 'Box', 'off');
% ax = gca;
% ax.YColor = 'k'; 
% ax.XColor = 'k'; % 保持刻度数字
% ax.TickDir = 'out';
% 




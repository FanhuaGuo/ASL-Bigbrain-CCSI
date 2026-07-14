clc
clear all
close all
set(0,'defaultfigurecolor',[1 1 1]);% set figure background color
if ~contains(path, '/Users/guofanhua/Desktop/gfh/tools/MatlabTools/BrainSpace-latest/matlab')
    addpath(genpath('/Users/guofanhua/Desktop/gfh/tools/MatlabTools/BrainSpace-latest/matlab'));
end

%% ======================= step 1 for disposal all data

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
    fprintf('Gene set number(Oligo 1->6) without end-states overlap:\n');
    for i = 1:size(CT_categories3,1)
        for j = 1:numel(CT_categories3.gene_label{i})
            if sum(tmpGene(CT_categories3.gene_label{i}(j),[1 5])) >= 1 && sum(tmpGene(CT_categories3.gene_label{i}(j),[2 3 4 6])) >= 1
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



%% compare retest-CCSI map and Mito-CI/MRC map 
CorrType = 'Pearson';  % Pearson Spearman
RepeatTimes = 20;
OutDir = '../Data/group/verification';
clear MitoCorr
mask = CBFsnr_mean_bil>4;
nullType = '_variogram';  %  _variogram  _spin  _eigenstrapping

% 5 subjects
ifload = 1;
PickNum = 5;
if ifload
    load([OutDir '/CCSI_5_Mito.mat']);
    MitoCorr.CI_5r = CI_r;
    MitoCorr.CI_5p = CI_p;
    MitoCorr.MRC_5r = MRC_r;
    MitoCorr.MRC_5p = MRC_p;
else
    CI_r = [];
    CI_p = [];
    MRC_r = [];
    MRC_p = [];
    for i = 1:RepeatTimes
        Data_name = [OutDir '/CCSI_' num2str(PickNum) '_' num2str(i) '_bil'];
        origx = load([Data_name '.1D']);
        origData1_null = load([Data_name nullType '.1D']);

        y = MitoData(:,1);
        x = origx(mask);
        y = y(mask);
        Data1_null = origData1_null(mask,:);
        c = corr(y,x,'Type',CorrType);
        cn = corr(Data1_null,y,'Type',CorrType);
        p = sum(c<cn)/20000;
        CI_r = [CI_r; c];
        CI_p = [CI_p; p];

        y = MitoData(:,5);
        x = origx(mask);
        y = y(mask);
        Data1_null = origData1_null(mask,:);
        c = corr(y,x,'Type',CorrType);
        cn = corr(Data1_null,y,'Type',CorrType);
        p = sum(c<cn)/20000;
        MRC_r = [MRC_r; c];
        MRC_p = [MRC_p; p];
    end
    save([OutDir '/CCSI_5_Mito.mat'],'CI_r','CI_p','MRC_r','MRC_p');
    MitoCorr.CI_5r = CI_r;
    MitoCorr.CI_5p = CI_p;
    MitoCorr.MRC_5r = MRC_r;
    MitoCorr.MRC_5p = MRC_p;
end



% 10 subjects
ifload = 1;
PickNum = 10;
if ifload
    load([OutDir '/CCSI_10_Mito.mat']);
    MitoCorr.CI_10r = CI_r;
    MitoCorr.CI_10p = CI_p;
    MitoCorr.MRC_10r = MRC_r;
    MitoCorr.MRC_10p = MRC_p;
else
    CI_r = [];
    CI_p = [];
    MRC_r = [];
    MRC_p = [];
    for i = 1:RepeatTimes
        Data_name = [OutDir '/CCSI_' num2str(PickNum) '_' num2str(i) '_bil'];
        origx = load([Data_name '.1D']);
        origData1_null = load([Data_name nullType '.1D']);

        y = MitoData(:,1);
        x = origx(mask);
        y = y(mask);
        Data1_null = origData1_null(mask,:);
        c = corr(y,x,'Type',CorrType);
        cn = corr(Data1_null,y,'Type',CorrType);
        p = sum(c<cn)/20000;
        CI_r = [CI_r; c];
        CI_p = [CI_p; p];

        y = MitoData(:,5);
        x = origx(mask);
        y = y(mask);
        Data1_null = origData1_null(mask,:);
        c = corr(y,x,'Type',CorrType);
        cn = corr(Data1_null,y,'Type',CorrType);
        p = sum(c<cn)/20000;
        MRC_r = [MRC_r; c];
        MRC_p = [MRC_p; p];
    end
    save([OutDir '/CCSI_10_Mito.mat'],'CI_r','CI_p','MRC_r','MRC_p');
    MitoCorr.CI_10r = CI_r;
    MitoCorr.CI_10p = CI_p;
    MitoCorr.MRC_10r = MRC_r;
    MitoCorr.MRC_10p = MRC_p;
end




% 15 subjects
ifload = 1;
PickNum = 15;
if ifload
    load([OutDir '/CCSI_15_Mito.mat']);
    MitoCorr.CI_15r = CI_r;
    MitoCorr.CI_15p = CI_p;
    MitoCorr.MRC_15r = MRC_r;
    MitoCorr.MRC_15p = MRC_p;
else
    CI_r = [];
    CI_p = [];
    MRC_r = [];
    MRC_p = [];
    for i = 1:RepeatTimes
        Data_name = [OutDir '/CCSI_' num2str(PickNum) '_' num2str(i) '_bil'];
        origx = load([Data_name '.1D']);
        origData1_null = load([Data_name nullType '.1D']);

        y = MitoData(:,1);
        x = origx(mask);
        y = y(mask);
        Data1_null = origData1_null(mask,:);
        c = corr(y,x,'Type',CorrType);
        cn = corr(Data1_null,y,'Type',CorrType);
        p = sum(c<cn)/20000;
        CI_r = [CI_r; c];
        CI_p = [CI_p; p];

        y = MitoData(:,5);
        x = origx(mask);
        y = y(mask);
        Data1_null = origData1_null(mask,:);
        c = corr(y,x,'Type',CorrType);
        cn = corr(Data1_null,y,'Type',CorrType);
        p = sum(c<cn)/20000;
        MRC_r = [MRC_r; c];
        MRC_p = [MRC_p; p];
    end
    save([OutDir '/CCSI_15_Mito.mat'],'CI_r','CI_p','MRC_r','MRC_p');
    MitoCorr.CI_15r = CI_r;
    MitoCorr.CI_15p = CI_p;
    MitoCorr.MRC_15r = MRC_r;
    MitoCorr.MRC_15p = MRC_p;
end



% 20 subjects
ifload = 1;
PickNum = 20;
if ifload
    load([OutDir '/CCSI_20_Mito.mat']);
    MitoCorr.CI_20r = CI_r;
    MitoCorr.CI_20p = CI_p;
    MitoCorr.MRC_20r = MRC_r;
    MitoCorr.MRC_20p = MRC_p;
else
    CI_r = [];
    CI_p = [];
    MRC_r = [];
    MRC_p = [];
    for i = 1:RepeatTimes
        Data_name = [OutDir '/CCSI_' num2str(PickNum) '_' num2str(i) '_bil'];
        origx = load([Data_name '.1D']);
        origData1_null = load([Data_name nullType '.1D']);

        y = MitoData(:,1);
        x = origx(mask);
        y = y(mask);
        Data1_null = origData1_null(mask,:);
        c = corr(y,x,'Type',CorrType);
        cn = corr(Data1_null,y,'Type',CorrType);
        p = sum(c<cn)/20000;
        CI_r = [CI_r; c];
        CI_p = [CI_p; p];

        y = MitoData(:,5);
        x = origx(mask);
        y = y(mask);
        Data1_null = origData1_null(mask,:);
        c = corr(y,x,'Type',CorrType);
        cn = corr(Data1_null,y,'Type',CorrType);
        p = sum(c<cn)/20000;
        MRC_r = [MRC_r; c];
        MRC_p = [MRC_p; p];
    end
    save([OutDir '/CCSI_20_Mito.mat'],'CI_r','CI_p','MRC_r','MRC_p');
    MitoCorr.CI_20r = CI_r;
    MitoCorr.CI_20p = CI_p;
    MitoCorr.MRC_20r = MRC_r;
    MitoCorr.MRC_20p = MRC_p;
end


%% compute CCSI CellType1-GSEA by median correlation's r-fisherZ
% common parameters
cate = CT_categories3.gene_label;
cateLabel = CT_categories3.Cell;
CorrType = 'Pearson';  % Pearson Spearman
GeneMask = GeneMask_SNR4;
Sample_Mask = Sample_Mask_SNR4;
nullType = '_variogram';  %  _variogram  _spin  _eigenstrapping
clear OligoEnrich

% 5 subjects
ifload = 1;
PickNum = 5;
if ifload
    load([OutDir '/CCSI_5_Oligo.mat']);
    OligoEnrich.Oligo1_5score = Oligo1_score;
    OligoEnrich.Oligo1_5p = Oligo1_p;
    OligoEnrich.Oligo5_5score = Oligo5_score;
    OligoEnrich.Oligo5_5p = Oligo5_p;
else
    Oligo1_score = [];
    Oligo1_p = [];
    Oligo5_score = [];
    Oligo5_p = [];
    for i = 1:RepeatTimes
        Data_name = [OutDir '/CCSI_' num2str(PickNum) '_' num2str(i) '_uni'];
        origx = load([Data_name '.1D']);
        origData1_null = load([Data_name nullType '.1D']);

        y = GeneExpression_Orig(cate{1},GeneMask)';
        x = origx(Sample_Mask);
        Data1_null = origData1_null(Sample_Mask,:);
        c = corr(y,x,'Type',CorrType);
        cn = corr(Data1_null,y,'Type',CorrType);
        cate_rZ = median(c);
        cate_rnZ = median(cn,2);
        p = (1+sum(abs(cate_rZ)<abs(cate_rnZ)))/20001;
        Oligo1_score = [Oligo1_score; cate_rZ];
        Oligo1_p = [Oligo1_p; p];

        y = GeneExpression_Orig(cate{5},GeneMask)';
        x = origx(Sample_Mask);
        Data1_null = origData1_null(Sample_Mask,:);
        c = corr(y,x,'Type',CorrType);
        cn = corr(Data1_null,y,'Type',CorrType);
        cate_rZ = median(c);
        cate_rnZ = median(cn,2);
        p = (1+sum(abs(cate_rZ)<abs(cate_rnZ)))/20001;
        Oligo5_score = [Oligo5_score; cate_rZ];
        Oligo5_p = [Oligo5_p; p];
    end
    save([OutDir '/CCSI_5_Oligo.mat'],'Oligo1_score','Oligo1_p','Oligo5_score','Oligo5_p');
    OligoEnrich.Oligo1_5score = Oligo1_score;
    OligoEnrich.Oligo1_5p = Oligo1_p;
    OligoEnrich.Oligo5_5score = Oligo5_score;
    OligoEnrich.Oligo5_5p = Oligo5_p;
end




% 10 subjects
ifload = 1;
PickNum = 10;
if ifload
    load([OutDir '/CCSI_10_Oligo.mat']);
    OligoEnrich.Oligo1_10score = Oligo1_score;
    OligoEnrich.Oligo1_10p = Oligo1_p;
    OligoEnrich.Oligo5_10score = Oligo5_score;
    OligoEnrich.Oligo5_10p = Oligo5_p;
else
    Oligo1_score = [];
    Oligo1_p = [];
    Oligo5_score = [];
    Oligo5_p = [];
    for i = 1:RepeatTimes
        Data_name = [OutDir '/CCSI_' num2str(PickNum) '_' num2str(i) '_uni'];
        origx = load([Data_name '.1D']);
        origData1_null = load([Data_name nullType '.1D']);

        y = GeneExpression_Orig(cate{1},GeneMask)';
        x = origx(Sample_Mask);
        Data1_null = origData1_null(Sample_Mask,:);
        c = corr(y,x,'Type',CorrType);
        cn = corr(Data1_null,y,'Type',CorrType);
        cate_rZ = median(c);
        cate_rnZ = median(cn,2);
        p = (1+sum(abs(cate_rZ)<abs(cate_rnZ)))/20001;
        Oligo1_score = [Oligo1_score; cate_rZ];
        Oligo1_p = [Oligo1_p; p];

        y = GeneExpression_Orig(cate{5},GeneMask)';
        x = origx(Sample_Mask);
        Data1_null = origData1_null(Sample_Mask,:);
        c = corr(y,x,'Type',CorrType);
        cn = corr(Data1_null,y,'Type',CorrType);
        cate_rZ = median(c);
        cate_rnZ = median(cn,2);
        p = (1+sum(abs(cate_rZ)<abs(cate_rnZ)))/20001;
        Oligo5_score = [Oligo5_score; cate_rZ];
        Oligo5_p = [Oligo5_p; p];
    end
    save([OutDir '/CCSI_10_Oligo.mat'],'Oligo1_score','Oligo1_p','Oligo5_score','Oligo5_p');
    OligoEnrich.Oligo1_10score = Oligo1_score;
    OligoEnrich.Oligo1_10p = Oligo1_p;
    OligoEnrich.Oligo5_10score = Oligo5_score;
    OligoEnrich.Oligo5_10p = Oligo5_p;
end



% 15 subjects
ifload = 1;
PickNum = 15;
if ifload
    load([OutDir '/CCSI_15_Oligo.mat']);
    OligoEnrich.Oligo1_15score = Oligo1_score;
    OligoEnrich.Oligo1_15p = Oligo1_p;
    OligoEnrich.Oligo5_15score = Oligo5_score;
    OligoEnrich.Oligo5_15p = Oligo5_p;
else
    Oligo1_score = [];
    Oligo1_p = [];
    Oligo5_score = [];
    Oligo5_p = [];
    for i = 1:RepeatTimes
        Data_name = [OutDir '/CCSI_' num2str(PickNum) '_' num2str(i) '_uni'];
        origx = load([Data_name '.1D']);
        origData1_null = load([Data_name nullType '.1D']);

        y = GeneExpression_Orig(cate{1},GeneMask)';
        x = origx(Sample_Mask);
        Data1_null = origData1_null(Sample_Mask,:);
        c = corr(y,x,'Type',CorrType);
        cn = corr(Data1_null,y,'Type',CorrType);
        cate_rZ = median(c);
        cate_rnZ = median(cn,2);
        p = (1+sum(abs(cate_rZ)<abs(cate_rnZ)))/20001;
        Oligo1_score = [Oligo1_score; cate_rZ];
        Oligo1_p = [Oligo1_p; p];

        y = GeneExpression_Orig(cate{5},GeneMask)';
        x = origx(Sample_Mask);
        Data1_null = origData1_null(Sample_Mask,:);
        c = corr(y,x,'Type',CorrType);
        cn = corr(Data1_null,y,'Type',CorrType);
        cate_rZ = median(c);
        cate_rnZ = median(cn,2);
        p = (1+sum(abs(cate_rZ)<abs(cate_rnZ)))/20001;
        Oligo5_score = [Oligo5_score; cate_rZ];
        Oligo5_p = [Oligo5_p; p];
    end
    save([OutDir '/CCSI_15_Oligo.mat'],'Oligo1_score','Oligo1_p','Oligo5_score','Oligo5_p');
    OligoEnrich.Oligo1_15score = Oligo1_score;
    OligoEnrich.Oligo1_15p = Oligo1_p;
    OligoEnrich.Oligo5_15score = Oligo5_score;
    OligoEnrich.Oligo5_15p = Oligo5_p;
end



% 20 subjects
ifload = 1;
PickNum = 20;
if ifload
    load([OutDir '/CCSI_20_Oligo.mat']);
    OligoEnrich.Oligo1_20score = Oligo1_score;
    OligoEnrich.Oligo1_20p = Oligo1_p;
    OligoEnrich.Oligo5_20score = Oligo5_score;
    OligoEnrich.Oligo5_20p = Oligo5_p;
else
    Oligo1_score = [];
    Oligo1_p = [];
    Oligo5_score = [];
    Oligo5_p = [];
    for i = 1:RepeatTimes
        Data_name = [OutDir '/CCSI_' num2str(PickNum) '_' num2str(i) '_uni'];
        origx = load([Data_name '.1D']);
        origData1_null = load([Data_name nullType '.1D']);

        y = GeneExpression_Orig(cate{1},GeneMask)';
        x = origx(Sample_Mask);
        Data1_null = origData1_null(Sample_Mask,:);
        c = corr(y,x,'Type',CorrType);
        cn = corr(Data1_null,y,'Type',CorrType);
        cate_rZ = median(c);
        cate_rnZ = median(cn,2);
        p = (1+sum(abs(cate_rZ)<abs(cate_rnZ)))/20001;
        Oligo1_score = [Oligo1_score; cate_rZ];
        Oligo1_p = [Oligo1_p; p];

        y = GeneExpression_Orig(cate{5},GeneMask)';
        x = origx(Sample_Mask);
        Data1_null = origData1_null(Sample_Mask,:);
        c = corr(y,x,'Type',CorrType);
        cn = corr(Data1_null,y,'Type',CorrType);
        cate_rZ = median(c);
        cate_rnZ = median(cn,2);
        p = (1+sum(abs(cate_rZ)<abs(cate_rnZ)))/20001;
        Oligo5_score = [Oligo5_score; cate_rZ];
        Oligo5_p = [Oligo5_p; p];
    end
    save([OutDir '/CCSI_20_Oligo.mat'],'Oligo1_score','Oligo1_p','Oligo5_score','Oligo5_p');
    OligoEnrich.Oligo1_20score = Oligo1_score;
    OligoEnrich.Oligo1_20p = Oligo1_p;
    OligoEnrich.Oligo5_20score = Oligo5_score;
    OligoEnrich.Oligo5_20p = Oligo5_p;
end


%% plotting mito-CI
% manual parameters
DataPlot = [[1*ones(RepeatTimes,1); 3*ones(RepeatTimes,1); 5*ones(RepeatTimes,1); 7*ones(RepeatTimes,1);]...
            [MitoCorr.CI_5r; MitoCorr.CI_10r; MitoCorr.CI_15r; MitoCorr.CI_20r]...
            [MitoCorr.CI_5p; MitoCorr.CI_10p; MitoCorr.CI_15p; MitoCorr.CI_20p]];
DataPlot(:,3) = -log10(DataPlot(:,3));
threshold = -log10([0.05 0.001]);

% color score
Plotscore = DataPlot(:,3);
Plotscore(Plotscore>threshold(2)) = threshold(2);
Plotscore(Plotscore<threshold(1)) = threshold(1);
Plotscore = 1+round(1000 * (Plotscore - threshold(1)) / (threshold(2) - threshold(1)));

% color gradient
colorGra = [ones(1,1000); 0.9:-0.9/999:0; 0.9:-0.9/999:0]';
colorGra = [0.6 0.6 0.6; colorGra];


% plotting
figure('Color',[1 1 1],'Position',[0 0 200 200],'Units','pixels');
hold on;
for i = 1:size(DataPlot)
    scatter(DataPlot(i,1), DataPlot(i,2), 50, colorGra(Plotscore(i),:), 'filled');
end
set(gca,'xLim',[0 8],'xTick',[1:2:7],'xTickLabel',{'R5','R10','R15','R20'}, 'FontName', 'Arial');
set(gca, 'ylim', [0.1 0.6], 'yTick', [0.1 0.6], 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';



%% plotting mito-MRC
% manual parameters
DataPlot = [[1*ones(RepeatTimes,1); 3*ones(RepeatTimes,1); 5*ones(RepeatTimes,1); 7*ones(RepeatTimes,1);]...
            [MitoCorr.MRC_5r; MitoCorr.MRC_10r; MitoCorr.MRC_15r; MitoCorr.MRC_20r]...
            [MitoCorr.MRC_5p; MitoCorr.MRC_10p; MitoCorr.MRC_15p; MitoCorr.MRC_20p]];
DataPlot(:,3) = -log10(DataPlot(:,3));
threshold = -log10([0.05 0.001]);

% color score
Plotscore = DataPlot(:,3);
Plotscore(Plotscore>threshold(2)) = threshold(2);
Plotscore(Plotscore<threshold(1)) = threshold(1);
Plotscore = 1+round(1000 * (Plotscore - threshold(1)) / (threshold(2) - threshold(1)));

% color gradient
colorGra = [ones(1,1000); 0.9:-0.9/999:0; 0.9:-0.9/999:0]';
colorGra = [0.6 0.6 0.6; colorGra];


% plotting
figure('Color',[1 1 1],'Position',[0 0 200 200],'Units','pixels');
hold on;
for i = 1:size(DataPlot)
    scatter(DataPlot(i,1), DataPlot(i,2), 50, colorGra(Plotscore(i),:), 'filled');
end
set(gca,'xLim',[0 8],'xTick',[1:2:7],'xTickLabel',{'R5','R10','R15','R20'}, 'FontName', 'Arial');
set(gca, 'ylim', [0.1 0.6], 'yTick', [0.1 0.6], 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';


        
        
%% plotting Oligo1
% manual parameters
DataPlot = [[1*ones(RepeatTimes,1); 3*ones(RepeatTimes,1); 5*ones(RepeatTimes,1); 7*ones(RepeatTimes,1);]...
            [OligoEnrich.Oligo1_5score; OligoEnrich.Oligo1_10score; OligoEnrich.Oligo1_15score; OligoEnrich.Oligo1_20score]...
            [OligoEnrich.Oligo1_5p; OligoEnrich.Oligo1_10p; OligoEnrich.Oligo1_15p; OligoEnrich.Oligo1_20p]];
DataPlot(:,3) = -log10(DataPlot(:,3));
threshold = -log10([0.05 0.0001]);

% color score
Plotscore = DataPlot(:,3);
Plotscore(Plotscore>threshold(2)) = threshold(2);
Plotscore(Plotscore<threshold(1)) = threshold(1);
Plotscore = 1+round(1000 * (Plotscore - threshold(1)) / (threshold(2) - threshold(1)));

% color gradient
colorGra = [ones(1,1000); 0.9:-0.9/999:0; 0.9:-0.9/999:0]';
colorGra = [0.6 0.6 0.6; colorGra];


% plotting
figure('Color',[1 1 1],'Position',[0 0 200 200],'Units','pixels');
hold on;
for i = 1:size(DataPlot)
    scatter(DataPlot(i,1), DataPlot(i,2), 50, colorGra(Plotscore(i),:), 'filled');
end
set(gca,'xLim',[0 8],'xTick',[1:2:7],'xTickLabel',{'R5','R10','R15','R20'}, 'FontName', 'Arial');
set(gca, 'ylim', [0.2 0.4], 'yTick', [0.2 0.4], 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';

        
%% plotting Oligo5
% manual parameters
DataPlot = [[1*ones(RepeatTimes,1); 3*ones(RepeatTimes,1); 5*ones(RepeatTimes,1); 7*ones(RepeatTimes,1);]...
            [OligoEnrich.Oligo5_5score; OligoEnrich.Oligo5_10score; OligoEnrich.Oligo5_15score; OligoEnrich.Oligo5_20score]...
            [OligoEnrich.Oligo5_5p; OligoEnrich.Oligo5_10p; OligoEnrich.Oligo5_15p; OligoEnrich.Oligo5_20p]];
DataPlot(:,3) = -log10(DataPlot(:,3));
threshold = -log10([0.05 0.0001]);

% color score
Plotscore = DataPlot(:,3);
Plotscore(Plotscore>threshold(2)) = threshold(2);
Plotscore(Plotscore<threshold(1)) = threshold(1);
Plotscore = 1+round(1000 * (Plotscore - threshold(1)) / (threshold(2) - threshold(1)));

% color gradient
colorGra = [ones(1,1000); 0.9:-0.9/999:0; 0.9:-0.9/999:0]';
colorGra = [0.6 0.6 0.6; colorGra];


% plotting
figure('Color',[1 1 1],'Position',[0 0 200 200],'Units','pixels');
hold on;
for i = 1:size(DataPlot)
    scatter(DataPlot(i,1), DataPlot(i,2), 50, colorGra(Plotscore(i),:), 'filled');
end
set(gca,'xLim',[0 8],'xTick',[1:2:7],'xTickLabel',{'R5','R10','R15','R20'}, 'FontName', 'Arial');
set(gca, 'ylim', [0.2 0.4], 'yTick', [0.2 0.4], 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';

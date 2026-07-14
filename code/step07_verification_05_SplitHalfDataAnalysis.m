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
OutDir = '../Data/group/SplitHalf';
clear MitoCorr MapCorr
mask = CBFsnr_mean_bil>4;
nullType = '_variogram';  %  _variogram  _spin  _eigenstrapping

% run
ifload = 1;
if ifload
    load([OutDir '/CCSI_Mito.mat']);
    MitoCorr.CI_r = CI_r;
    MitoCorr.CI_p = CI_p;
    MitoCorr.MRC_r = MRC_r;
    MitoCorr.MRC_p = MRC_p;
    load([OutDir '/CCSI_MapCorr.mat']);
    MapCorr.r = Map_r;
    MapCorr.p = Map_p;
else
    CI_r = [];
    CI_p = [];
    MRC_r = [];
    MRC_p = [];
    Map_r = [];
    Map_p = [];
    for i = 1:RepeatTimes
        Data_name = [OutDir '/CCSI1_' num2str(i) '_bil'];
        origx1 = load([Data_name '.1D']);
        origData1_null = load([Data_name nullType '.1D']);
        Data_name = [OutDir '/CCSI2_' num2str(i) '_bil'];
        origx2 = load([Data_name '.1D']);
        origData2_null = load([Data_name nullType '.1D']);

        % CI
        y = MitoData(:,1);
        x = origx1(mask);
        y = y(mask);
        Data1_null = origData1_null(mask,:);
        c1 = corr(y,x,'Type',CorrType);
        cn1 = corr(Data1_null,y,'Type',CorrType);
        p1 = sum(abs(c1)<abs(cn1))/20000;
        
        y = MitoData(:,1);
        x = origx2(mask);
        y = y(mask);
        Data2_null = origData2_null(mask,:);
        c2 = corr(y,x,'Type',CorrType);
        cn2 = corr(Data1_null,y,'Type',CorrType);
        p2 = sum(abs(c2)<abs(cn2))/20000;
        
        CI_r = [CI_r; c1 c2];
        CI_p = [CI_p; p1 p2];
        
        % MRC
        y = MitoData(:,5);
        x = origx1(mask);
        y = y(mask);
        Data1_null = origData1_null(mask,:);
        c1 = corr(y,x,'Type',CorrType);
        cn1 = corr(Data1_null,y,'Type',CorrType);
        p1 = sum(abs(c1)<abs(cn1))/20000;
        
        y = MitoData(:,5);
        x = origx2(mask);
        y = y(mask);
        Data2_null = origData2_null(mask,:);
        c2 = corr(y,x,'Type',CorrType);
        cn2 = corr(Data1_null,y,'Type',CorrType);
        p2 = sum(abs(c2)<abs(cn2))/20000;
        
        MRC_r = [MRC_r; c1 c2];
        MRC_p = [MRC_p; p1 p2];
        
        % Map corr
        y = origx2(mask);
        x = origx1(mask);
        Data1_null = origData1_null(mask,:);
        c = corr(y,x,'Type',CorrType);
        cn = corr(Data1_null,y,'Type',CorrType);
        p = sum(abs(c)<abs(cn))/20000;
        
        Map_r = [Map_r; c];
        Map_p = [Map_p; p];
    end
    save([OutDir '/CCSI_Mito.mat'],'CI_r','CI_p','MRC_r','MRC_p');
    MitoCorr.CI_r = CI_r;
    MitoCorr.CI_p = CI_p;
    MitoCorr.MRC_r = MRC_r;
    MitoCorr.MRC_p = MRC_p;
    save([OutDir '/CCSI_MapCorr.mat'],'Map_r','Map_p');
    MapCorr.r = Map_r;
    MapCorr.p = Map_p;
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

% run
ifload = 1;
if ifload
    load([OutDir '/CCSI_Oligo.mat']);
    OligoEnrich.Oligo1_score = Oligo1_score;
    OligoEnrich.Oligo1_p = Oligo1_p;
    OligoEnrich.Oligo5_score = Oligo5_score;
    OligoEnrich.Oligo5_p = Oligo5_p;
else
    Oligo1_score = [];
    Oligo1_p = [];
    Oligo5_score = [];
    Oligo5_p = [];
    for i = 1:RepeatTimes        
        Data_name = [OutDir '/CCSI1_' num2str(i) '_uni'];
        origx1 = load([Data_name '.1D']);
        origData1_null = load([Data_name nullType '.1D']);
        Data_name = [OutDir '/CCSI2_' num2str(i) '_uni'];
        origx2 = load([Data_name '.1D']);
        origData2_null = load([Data_name nullType '.1D']);

        % oligo1
        y = GeneExpression_Orig(cate{1},GeneMask)';
        x = origx1(Sample_Mask);
        Data1_null = origData1_null(Sample_Mask,:);
        c = corr(y,x,'Type',CorrType);
        cn = corr(Data1_null,y,'Type',CorrType);
        cate_rZ1 = median(c);
        cate_rnZ1 = median(cn,2);
        p1 = (1+sum(abs(cate_rZ1)<abs(cate_rnZ1)))/20001;
        
        y = GeneExpression_Orig(cate{1},GeneMask)';
        x = origx2(Sample_Mask);
        Data2_null = origData2_null(Sample_Mask,:);
        c = corr(y,x,'Type',CorrType);
        cn = corr(Data2_null,y,'Type',CorrType);
        cate_rZ2 = median(c);
        cate_rnZ2 = median(cn,2);
        p2 = (1+sum(abs(cate_rZ2)<abs(cate_rnZ2)))/20001;
        
        Oligo1_score = [Oligo1_score; cate_rZ1 cate_rZ2];
        Oligo1_p = [Oligo1_p; p1 p2];
        
        % oligo5
        y = GeneExpression_Orig(cate{5},GeneMask)';
        x = origx1(Sample_Mask);
        Data1_null = origData1_null(Sample_Mask,:);
        c = corr(y,x,'Type',CorrType);
        cn = corr(Data1_null,y,'Type',CorrType);
        cate_rZ1 = median(c);
        cate_rnZ1 = median(cn,2);
        p1 = (1+sum(abs(cate_rZ1)<abs(cate_rnZ1)))/20001;
        
        y = GeneExpression_Orig(cate{5},GeneMask)';
        x = origx2(Sample_Mask);
        Data2_null = origData2_null(Sample_Mask,:);
        c = corr(y,x,'Type',CorrType);
        cn = corr(Data2_null,y,'Type',CorrType);
        cate_rZ2 = median(c);
        cate_rnZ2 = median(cn,2);
        p2 = (1+sum(abs(cate_rZ2)<abs(cate_rnZ2)))/20001;
        
        Oligo5_score = [Oligo5_score; cate_rZ1 cate_rZ2];
        Oligo5_p = [Oligo5_p; p1 p2];
    end
    save([OutDir '/CCSI_Oligo.mat'],'Oligo1_score','Oligo1_p','Oligo5_score','Oligo5_p');
    OligoEnrich.Oligo1_score = Oligo1_score;
    OligoEnrich.Oligo1_p = Oligo1_p;
    OligoEnrich.Oligo5_score = Oligo5_score;
    OligoEnrich.Oligo5_p = Oligo5_p;
end



%% plotting a correlation scatter
Which_time = 1;
x = load([OutDir '/CCSI1_' num2str(Which_time) '_bil.1D']);
y = load([OutDir '/CCSI2_' num2str(Which_time) '_bil.1D']);
r = MapCorr.r(Which_time);
p = MapCorr.p(Which_time);
fprintf(['r = ' num2str(r) ', variogram-p = ' num2str(p) '\n']);

% plot
figure('Color',[1 1 1],'Position',[0 0 300 300],'Units','pixels');
hold on;
subfun_plot_scatter_regression( x , y , [-0.5 1.2] , [-0.5 1.2]);


% write to SourceData Fig7h
outSD = [x y];


%% plotting analysis
% manual parameters
threshold = -log10([0.05 0.0001]);
markerSize = 80;

% write to SourceData Fig7i
outSD = [];

% color gradient
colorGra = [ones(1,1000); 0.9:-0.9/999:0; 0.9:-0.9/999:0]';
colorGra = [0.6 0.6 0.6; colorGra];


% plotting
figure('Color',[1 1 1],'Position',[0 0 600 300],'Units','pixels');
hold on;

% spatial corr
r = MapCorr.r;
% r = (r-0.8)*2+0.6;
r = r-0.2;
tmp = MapCorr.p;
tmp(tmp==0) = 1/20001;
Plotscore = -log10(tmp); % color score
Plotscore(Plotscore>threshold(2)) = threshold(2);
Plotscore(Plotscore<threshold(1)) = threshold(1);
Plotscore = 1+round(1000 * (Plotscore - threshold(1)) / (threshold(2) - threshold(1)));
scatter(-1*ones(RepeatTimes,1), r, markerSize, colorGra(Plotscore,:), '^', 'filled'); % plot
outSD = [outSD MapCorr.r tmp];

% CI
x = [3 4];
r = MitoCorr.CI_r;
tmp = MitoCorr.CI_p;
tmp(tmp==0) = 1/20001;
Plotscore = -log10(tmp); % color score
Plotscore(Plotscore>threshold(2)) = threshold(2);
Plotscore(Plotscore<threshold(1)) = threshold(1);
Plotscore = 1+round(1000 * (Plotscore - threshold(1)) / (threshold(2) - threshold(1)));
scatter(x(1)*ones(RepeatTimes,1), r(:,1), markerSize, colorGra(Plotscore(:,1),:), 'o', 'filled'); % plot
scatter(x(2)*ones(RepeatTimes,1), r(:,2), markerSize, colorGra(Plotscore(:,2),:), 's', 'filled'); % plot
outSD = [outSD r(:,1) tmp(:,1) r(:,2) tmp(:,2)];

% MRC
x = [6 7];
r = MitoCorr.MRC_r;
tmp = MitoCorr.MRC_p;
tmp(tmp==0) = 1/20001;
Plotscore = -log10(tmp); % color score
Plotscore(Plotscore>threshold(2)) = threshold(2);
Plotscore(Plotscore<threshold(1)) = threshold(1);
Plotscore = 1+round(1000 * (Plotscore - threshold(1)) / (threshold(2) - threshold(1)));
scatter(x(1)*ones(RepeatTimes,1), r(:,1), markerSize, colorGra(Plotscore(:,1),:), 'o', 'filled'); % plot
scatter(x(2)*ones(RepeatTimes,1), r(:,2), markerSize, colorGra(Plotscore(:,2),:), 's', 'filled'); % plot
outSD = [outSD r(:,1) tmp(:,1) r(:,2) tmp(:,2)];

% oligo1
x = [9 10];
r = OligoEnrich.Oligo1_score;
tmp = OligoEnrich.Oligo1_p;
tmp(tmp==0) = 1/20001;
Plotscore = -log10(tmp); % color score
Plotscore(Plotscore>threshold(2)) = threshold(2);
Plotscore(Plotscore<threshold(1)) = threshold(1);
Plotscore = 1+round(1000 * (Plotscore - threshold(1)) / (threshold(2) - threshold(1)));
scatter(x(1)*ones(RepeatTimes,1), r(:,1), markerSize, colorGra(Plotscore(:,1),:), 'o', 'filled'); % plot
scatter(x(2)*ones(RepeatTimes,1), r(:,2), markerSize, colorGra(Plotscore(:,2),:), 's', 'filled'); % plot
outSD = [outSD r(:,1) tmp(:,1) r(:,2) tmp(:,2)];

% oligo5
x = [12 13];
r = OligoEnrich.Oligo5_score;
tmp = OligoEnrich.Oligo5_p;
tmp(tmp==0) = 1/20001;
Plotscore = -log10(tmp); % color score
Plotscore(Plotscore>threshold(2)) = threshold(2);
Plotscore(Plotscore<threshold(1)) = threshold(1);
Plotscore = 1+round(1000 * (Plotscore - threshold(1)) / (threshold(2) - threshold(1)));
scatter(x(1)*ones(RepeatTimes,1), r(:,1), markerSize, colorGra(Plotscore(:,1),:), 'o', 'filled'); % plot
scatter(x(2)*ones(RepeatTimes,1), r(:,2), markerSize, colorGra(Plotscore(:,2),:), 's', 'filled'); % plot
outSD = [outSD r(:,1) tmp(:,1) r(:,2) tmp(:,2)];

% set
set(gca,'xLim',[-2 14],'xTick',[-1 3.5 6.5 9.5 12.5],'xTickLabel',{'spatial corr','CI','MRC','oligo1','oligo5'}, 'FontName', 'Arial');
% set(gca, 'ylim', [0.2 1], 'yTick', [0.2 0.5 0.6 1], 'yTickLabel', [0.2 0.5 0.8 1], 'FontName', 'Arial');
set(gca, 'ylim', [0.2 0.8], 'yTick', [0.2 0.5 0.6 0.8], 'yTickLabel', [0.2 0.5 0.8 1], 'FontName', 'Arial');
set(gca, 'Box', 'off');
ax = gca;
ax.YColor = 'k'; 
ax.XColor = 'k'; % 保持刻度数字
ax.TickDir = 'out';




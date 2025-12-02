clc
clear all
close all
set(0,'defaultfigurecolor',[1 1 1]);% set figure background color
addpath(genpath('/Users/guofanhua/Desktop/gfh/tools/MatlabTools/BrainSpace-latest/matlab'))

%% ======================= compare regional CCSI

%% input data
step01_CBFandCCSImaps_00_init;
step01_CBFandCCSImaps_00_ComputeBasicResults


%% Mito maps info
info = {'CI' 'CII' 'CIV' 'MitoD' 'MRC' 'TRC'}';
for i = 1:6
    fprintf([info{i} ' Mean: ' num2str(mean(MitoData(:,i))) '  Std: ' num2str(std(MitoData(:,i))) '\n'])
end

%% set
CorrType = 'Pearson';  % Pearson Spearman


%% compare CBF-score/CCSI map and Mito-CI/CII/CIV map (Fig 2bc)
figure('Color',[1 1 1],'Position',[0 0 400 600],'Units','pixels');

% CBF-score
Data_name = 'CBF_PC1_Val_bil';  % 'CBF_PC1_Val_bil', 'CBF_Mean_PScs_bil'
origx = load([Group_dir Data_name '.1D']);
origData1_null = load([Group_dir Data_name '_null.1D']);
mask = CBFsnr_mean_bil>0;

r = [];
spinp = [];
for i = 1:3
    y = MitoData(:,i);
    x = origx(mask);
    y = y(mask);
    Data1_null = origData1_null(mask,:);

    c = corr(y,x,'Type',CorrType);
    cn = corr(Data1_null,y,'Type',CorrType);
    p = sum(c<cn)/20000;
    r = [r; c];
    spinp = [spinp; p];

    subplot(3,2,i*2-1); hold on;
    subfun_plot_scatter_regression( x , y , [-10 8] , [floor(min(y*10))/10 ceil(max(y*10))/10]);
end
r_Mito_complex = r;
p_Mito_complex = spinp;


% CCSI
Data_name = 'CBF_Mean_PScs_bil';  % 'CBF_PC1_Val_bil', 'CBF_Mean_PScs_bil'
origx = load([Group_dir Data_name '.1D']);
origData1_null = load([Group_dir Data_name '_null.1D']);
mask = CBFsnr_mean_bil>4;

r = [];
spinp = [];
for i = 1:3
    y = MitoData(:,i);
    x = origx(mask);
    y = y(mask);
    Data1_null = origData1_null(mask,:);

    c = corr(y,x,'Type',CorrType);
    cn = corr(Data1_null,y,'Type',CorrType);
    p = sum(c<cn)/20000;
    r = [r; c];
    spinp = [spinp; p];

    subplot(3,2,i*2); hold on;
    subfun_plot_scatter_regression( x , y , [0 1.2] , [floor(min(y*10))/10 ceil(max(y*10))/10]);
end
r_Mito_complex = [r_Mito_complex r];
p_Mito_complex = [p_Mito_complex spinp];

p_Mito_complex_FDR = fdr_BH(p_Mito_complex(:,1),0.05)'
p_Mito_complex_FDR = fdr_BH(p_Mito_complex(:,2),0.05)'





%% compare CBF-score/CCSI map and Mito-TRC/MRC map (Fig 2d and Fig S4)
figure('Color',[1 1 1],'Position',[0 0 400 600],'Units','pixels');

% CBF-score
Data_name = 'CBF_PC2_Val_bil';  % 'CBF_PC1_Val_bil', 'CBF_Mean_PScs_bil'
origx = load([Group_dir Data_name '.1D']);
origData1_null = load([Group_dir Data_name '_null.1D']);
mask = CBFsnr_mean_bil>0;

r = [];
spinp = [];
for i = 1:3
    y = MitoData(:,i+3);
    x = origx(mask);
    y = y(mask);
    Data1_null = origData1_null(mask,:);

    c = corr(y,x,'Type',CorrType);
    cn = corr(Data1_null,y,'Type',CorrType);
    p = sum(c<cn)/20000;
    r = [r; c];
    spinp = [spinp; p];

    subplot(3,2,i*2-1); hold on;
    subfun_plot_scatter_regression( x , y , [-10 8] , [floor(min(y*10))/10 ceil(max(y*10))/10]);
end
r_Mito_RC = r;
p_Mito_RC = spinp;


% CCSI
Data_name = 'CBF_Mean_PScs_bil';  % 'CBF_PC1_Val_bil', 'CBF_Mean_PScs_bil'
origx = load([Group_dir Data_name '.1D']);
origData1_null = load([Group_dir Data_name '_null.1D']);
mask = CBFsnr_mean_bil>4;

r = [];
spinp = [];
for i = 1:3
    y = MitoData(:,i+3);
    x = origx(mask);
    y = y(mask);
    Data1_null = origData1_null(mask,:);

    c = corr(y,x,'Type',CorrType);
    cn = corr(Data1_null,y,'Type',CorrType);
    p = sum(c<cn)/20000;
    r = [r; c];
    spinp = [spinp; p];

    subplot(3,2,i*2); hold on;
    subfun_plot_scatter_regression( x , y , [0 1.2] , [floor(min(y*10))/10 ceil(max(y*10))/10]);
end
r_Mito_RC = [r_Mito_RC r];
p_Mito_RC = [p_Mito_RC spinp];


p_Mito_RC_FDR = fdr_BH(p_Mito_RC(:,1),0.05)'
p_Mito_RC_FDR = fdr_BH(p_Mito_RC(:,2),0.05)'






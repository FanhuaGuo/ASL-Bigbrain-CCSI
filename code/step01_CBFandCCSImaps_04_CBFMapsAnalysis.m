clc
clear all
close all
set(0,'defaultfigurecolor',[1 1 1]);% set figure background color
addpath(genpath('/Users/guofanhua/Desktop/gfh/tools/MatlabTools/BrainSpace-latest/matlab'))

%% ======================= compare regional CCSI

%% input data
step01_CBFandCCSImaps_00_init;
step01_CBFandCCSImaps_00_ComputeBasicResults


%% compare CBF-score map and previous PET-CBF map (Fig S2a)
x = all_other_maps.maps{54};
Data_name = 'CBF_PC1_Val_bil';  % 'CBF_PC1_Val_bil', 'CBF_Mean_PScs_bil'
y = load([Group_dir Data_name '.1D']);
Data1_null = load([Group_dir Data_name '_null.1D']);
CorrType = 'Pearson';  % Pearson Spearman
ifdensity = 1;
ifisopotential = 1;
ifregressline = 0;
valueorrank = 2;  % 1:value 2:rank
ifDeOutlier = 0;

mask = CBFsnr_mean_bil>0 & abs(x-mean(x))<3*std(x);
x = x(mask);
y = y(mask);
Data1_null = Data1_null(mask,:);

c = corr(y,x,'Type',CorrType)
cn = corr(Data1_null,x,'Type',CorrType);
p = sum(abs(c)<abs(cn))/20000

if ifDeOutlier
    mask = subfun_FeatureSpace_DeOutlier(x);
    y = y(mask>0);
    x = x(mask>0);
end

figure('Color',[1 1 1],'Position',[0 0 300 300],'Units','pixels');
hold on;
subfun_plot_scatter_density_isopotential( x , y , valueorrank , ifdensity , ifisopotential , ifregressline )



%% compare CBF-score map and previous ASL-CBF map (Fig S2b)
x = all_other_maps.maps{64};
Data_name = 'CBF_PC1_Val_bil';  % 'CBF_PC1_Val_bil', 'CBF_Mean_PScs_bil'
y = load([Group_dir Data_name '.1D']);
Data1_null = load([Group_dir Data_name '_null.1D']);
CorrType = 'Pearson';  % Pearson Spearman
ifdensity = 1;
ifisopotential = 1;
ifregressline = 0;
valueorrank = 2;  % 1:value 2:rank
ifDeOutlier = 0;

mask = CBFsnr_mean_bil>0 & abs(x-mean(x))<3*std(x);
x = x(mask);
y = y(mask);
Data1_null = Data1_null(mask,:);

c = corr(y,x,'Type',CorrType)
cn = corr(Data1_null,x,'Type',CorrType);
p = sum(abs(c)<abs(cn))/20000

if ifDeOutlier
    mask = subfun_FeatureSpace_DeOutlier(x);
    y = y(mask>0);
    x = x(mask>0);
end

figure('Color',[1 1 1],'Position',[0 0 300 300],'Units','pixels');
hold on;
subfun_plot_scatter_density_isopotential( x , y , valueorrank , ifdensity , ifisopotential , ifregressline )




%% compare CBF-score map and previous PET-CMRO2 map (Fig S2c)
x = all_other_maps.maps{56};
Data_name = 'CBF_PC1_Val_bil';  % 'CBF_PC1_Val_bil', 'CBF_Mean_PScs_bil'
y = load([Group_dir Data_name '.1D']);
Data1_null = load([Group_dir Data_name '_null.1D']);
CorrType = 'Pearson';  % Pearson Spearman
ifdensity = 1;
ifisopotential = 1;
ifregressline = 0;
valueorrank = 2;  % 1:value 2:rank
ifDeOutlier = 0;

mask = CBFsnr_mean_bil>0 & abs(x-mean(x))<3*std(x);
x = x(mask);
y = y(mask);
Data1_null = Data1_null(mask,:);

c = corr(y,x,'Type',CorrType)
cn = corr(Data1_null,x,'Type',CorrType);
p = sum(abs(c)<abs(cn))/20000

if ifDeOutlier
    mask = subfun_FeatureSpace_DeOutlier(x);
    y = y(mask>0);
    x = x(mask>0);
end

figure('Color',[1 1 1],'Position',[0 0 300 300],'Units','pixels');
hold on;
subfun_plot_scatter_density_isopotential( x , y , valueorrank , ifdensity , ifisopotential , ifregressline )




%% compare CBF-score map and previous PET-CMRglc map (Fig S2d)
x = all_other_maps.maps{57};
Data_name = 'CBF_PC1_Val_bil';  % 'CBF_PC1_Val_bil', 'CBF_Mean_PScs_bil'
y = load([Group_dir Data_name '.1D']);
Data1_null = load([Group_dir Data_name '_null.1D']);
CorrType = 'Pearson';  % Pearson Spearman
ifdensity = 1;
ifisopotential = 1;
ifregressline = 0;
valueorrank = 2;  % 1:value 2:rank
ifDeOutlier = 0;

mask = CBFsnr_mean_bil>0 & abs(x-mean(x))<3*std(x);
x = x(mask);
y = y(mask);
Data1_null = Data1_null(mask,:);

c = corr(y,x,'Type',CorrType)
cn = corr(Data1_null,x,'Type',CorrType);
p = sum(abs(c)<abs(cn))/20000

if ifDeOutlier
    mask = subfun_FeatureSpace_DeOutlier(x);
    y = y(mask>0);
    x = x(mask>0);
end

figure('Color',[1 1 1],'Position',[0 0 300 300],'Units','pixels');
hold on;
subfun_plot_scatter_density_isopotential( x , y , valueorrank , ifdensity , ifisopotential , ifregressline )





%% compare CBF-score map and previous PET-martinostat(HDAC) map (Fig S2e)
x = all_other_maps.maps{75};
Data_name = 'CBF_PC1_Val_bil';  % 'CBF_PC1_Val_bil', 'CBF_Mean_PScs_bil'
y = load([Group_dir Data_name '.1D']);
Data1_null = load([Group_dir Data_name '_null.1D']);
CorrType = 'Pearson';  % Pearson Spearman
ifdensity = 1;
ifisopotential = 1;
ifregressline = 0;
valueorrank = 2;  % 1:value 2:rank
ifDeOutlier = 0;

mask = CBFsnr_mean_bil>0 & abs(x-mean(x))<3*std(x);
x = x(mask);
y = y(mask);
Data1_null = Data1_null(mask,:);

c = corr(y,x,'Type',CorrType)
cn = corr(Data1_null,x,'Type',CorrType);
p = sum(abs(c)<abs(cn))/20000

if ifDeOutlier
    mask = subfun_FeatureSpace_DeOutlier(x);
    y = y(mask>0);
    x = x(mask>0);
end

figure('Color',[1 1 1],'Position',[0 0 300 300],'Units','pixels');
hold on;
subfun_plot_scatter_density_isopotential( x , y , valueorrank , ifdensity , ifisopotential , ifregressline )



%% compare CBF-score map and previous scaling-PNC map (Fig S2e)
x = all_other_maps.maps{60};
Data_name = 'CBF_PC1_Val_bil';  % 'CBF_PC1_Val_bil', 'CBF_Mean_PScs_bil'
y = load([Group_dir Data_name '.1D']);
Data1_null = load([Group_dir Data_name '_null.1D']);
CorrType = 'Pearson';  % Pearson Spearman
ifdensity = 1;
ifisopotential = 1;
ifregressline = 0;
valueorrank = 2;  % 1:value 2:rank
ifDeOutlier = 0;

mask = CBFsnr_mean_bil>0 & abs(x-mean(x))<3*std(x);
x = x(mask);
y = y(mask);
Data1_null = Data1_null(mask,:);

c = corr(y,x,'Type',CorrType)
cn = corr(Data1_null,x,'Type',CorrType);
p = sum(abs(c)<abs(cn))/20000

if ifDeOutlier
    mask = subfun_FeatureSpace_DeOutlier(x);
    y = y(mask>0);
    x = x(mask>0);
end

figure('Color',[1 1 1],'Position',[0 0 300 300],'Units','pixels');
hold on;
subfun_plot_scatter_density_isopotential( x , y , valueorrank , ifdensity , ifisopotential , ifregressline )




clc
clear all
close all
set(0,'defaultfigurecolor',[1 1 1]);% set figure background color
if ~contains(path, '/Users/guofanhua/Desktop/gfh/tools/MatlabTools/BrainSpace-latest/matlab')
    addpath(genpath('/Users/guofanhua/Desktop/gfh/tools/MatlabTools/BrainSpace-latest/matlab'));
end

%% ======================= cbf maps analysis

%% input data
step01_CBFandCCSImaps_00_init;
step01_CBFandCCSImaps_00_ComputeBasicResults
corr(CBF_ValPC1_bil,CBF_Valmean_bil)

%% ======================= compare left and right hemi (symmetric analysis)
%% compare average lh and rh (Fig 2b)
CBF_lh_mean = mean(Data_Val_bil(1:180,:),1)';
CBF_rh_mean = mean(Data_Val_bil(181:360,:),1)';

% test 
[~,p_hemi] = ttest(CBF_lh_mean,CBF_rh_mean)

% plot
figure('Color',[1 1 1],'Position',[0 0 200 300],'Units','pixels');
hold on;

% 画每个被试的连线
for i = 1:numel(CBF_lh_mean)
    plot([1 2], [CBF_lh_mean(i) CBF_rh_mean(i)], '-o', ...
        'Color', [0.7 0.7 0.7], 'MarkerSize', 5, 'Linewidth', 0.5);
end

% 画均值（加粗）
plot([1 2], [mean(CBF_lh_mean) mean(CBF_rh_mean)], '-o', ...
    'Color', 'k', 'LineWidth', 2, 'MarkerSize', 8, 'Linewidth', 1.5);

xlim([0.5 2.5])
xticks([1 2])
xticklabels({'Left','Right'})
ylim([30,90]);
yticks([30:30:90]);
ylabel('CBF (mL/100g/min)')

title('Paired CBF (Left vs Right)')
box off


% write to SourceData Fig2b
outSD = [CBF_lh_mean CBF_rh_mean];


%% spatial distribution compare (asymmetric) (Fig 2c)
Data_name = ['CBF_Mean_Val_L12_bil'];  % 'CBF_PC1_Val_bil', 'CBF_Mean_PScs_bil'
y = load([Group_dir Data_name '.1D']);
Data1_null = load([Group_dir Data_name '_variogram.1D']);
CorrType = 'Pearson';  % Pearson Spearman

% plot
figure('Color',[1 1 1],'Position',[0 0 300 300],'Units','pixels');
hold on;
subfun_plot_scatter_regression(y(1:end/2) , y(end/2+1:end) , [30 80] , [30 80]);

% variogram p
c = corr(y(1:end/2), y(end/2+1:end),'Type',CorrType,'Rows', 'pairwise')
cn = corr(Data1_null(1:end/2,:), y(end/2+1:end),'Type',CorrType,'Rows', 'pairwise');
p = sum(abs(c)<abs(cn))/20000


% write to SourceData Fig2c
outSD = [atlas_info.region(1:180)];
tmp = table(y(1:end/2), y(end/2+1:end),...
            'VariableNames',{'CBF lh', 'CBF rh'});
outSD = [outSD tmp];



%% ======================= correlation with previous cbf maps
%% compare CBF-score map and previous PET-CBF map (Fig S2a)
x = all_other_maps.maps{54};
Data_name = 'CBF_PC1_Val_L12_bil';  % 'CBF_PC1_Val_bil', 'CBF_Mean_PScs_bil'
y = load([Group_dir Data_name '.1D']);
Data1_null = load([Group_dir Data_name '_variogram.1D']);
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

c = corr(y,x,'Type',CorrType,'Rows', 'pairwise')
cn = corr(Data1_null,x,'Type',CorrType,'Rows', 'pairwise');
p = sum(abs(c)<abs(cn))/20000

if ifDeOutlier
    mask = subfun_FeatureSpace_DeOutlier(x);
    y = y(mask>0);
    x = x(mask>0);
end

figure('Color',[1 1 1],'Position',[0 0 300 300],'Units','pixels');
hold on;
% subfun_plot_scatter_density_isopotential( x , y , valueorrank , ifdensity , ifisopotential , ifregressline )
subfun_plot_scatter_regression(y , x , [-12 9] , [4000 7000]);


%% compare CBF-score map and previous ASL-CBF map (Fig S2b)
x = all_other_maps.maps{64};
Data_name = 'CBF_PC1_Val_L12_bil';  % 'CBF_PC1_Val_bil', 'CBF_Mean_PScs_bil'
y = load([Group_dir Data_name '.1D']);
Data1_null = load([Group_dir Data_name '_variogram.1D']);
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

c = corr(y,x,'Type',CorrType,'Rows', 'pairwise')
cn = corr(Data1_null,x,'Type',CorrType,'Rows', 'pairwise');
p = sum(abs(c)<abs(cn))/20000

if ifDeOutlier
    mask = subfun_FeatureSpace_DeOutlier(x);
    y = y(mask>0);
    x = x(mask>0);
end

figure('Color',[1 1 1],'Position',[0 0 300 300],'Units','pixels');
hold on;
% subfun_plot_scatter_density_isopotential( x , y , valueorrank , ifdensity , ifisopotential , ifregressline )
subfun_plot_scatter_regression(y, x, [-12 9] , [40 90]);




%% compare CBF-score map and previous PET-CMRO2 map (Fig S2c)
x = all_other_maps.maps{56};
Data_name = 'CBF_PC1_Val_L12_bil';  % 'CBF_PC1_Val_bil', 'CBF_Mean_PScs_bil'
y = load([Group_dir Data_name '.1D']);
Data1_null = load([Group_dir Data_name '_variogram.1D']);
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

c = corr(y,x,'Type',CorrType,'Rows', 'pairwise')
cn = corr(Data1_null,x,'Type',CorrType,'Rows', 'pairwise');
p = sum(abs(c)<abs(cn))/20000

if ifDeOutlier
    mask = subfun_FeatureSpace_DeOutlier(x);
    y = y(mask>0);
    x = x(mask>0);
end

figure('Color',[1 1 1],'Position',[0 0 300 300],'Units','pixels');
hold on;
% subfun_plot_scatter_density_isopotential( x , y , valueorrank , ifdensity , ifisopotential , ifregressline )
subfun_plot_scatter_regression(y, x, [-12 9] , [4000 6000]);




%% compare CBF-score map and previous PET-CMRglc map (Fig S2d)
x = all_other_maps.maps{57};
Data_name = 'CBF_PC1_Val_L12_bil';  % 'CBF_PC1_Val_bil', 'CBF_Mean_PScs_bil'
y = load([Group_dir Data_name '.1D']);
Data1_null = load([Group_dir Data_name '_variogram.1D']);
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

c = corr(y,x,'Type',CorrType,'Rows', 'pairwise')
cn = corr(Data1_null,x,'Type',CorrType,'Rows', 'pairwise');
p = sum(abs(c)<abs(cn))/20000

if ifDeOutlier
    mask = subfun_FeatureSpace_DeOutlier(x);
    y = y(mask>0);
    x = x(mask>0);
end

figure('Color',[1 1 1],'Position',[0 0 300 300],'Units','pixels');
hold on;
% subfun_plot_scatter_density_isopotential( x , y , valueorrank , ifdensity , ifisopotential , ifregressline )
subfun_plot_scatter_regression(y, x, [-12 9] , [4000 8000]);





%% compare CBF-score map and previous PET-martinostat(HDAC) map (Fig S2e)
x = all_other_maps.maps{75};
Data_name = 'CBF_PC1_Val_L12_bil';  % 'CBF_PC1_Val_bil', 'CBF_Mean_PScs_bil'
y = load([Group_dir Data_name '.1D']);
Data1_null = load([Group_dir Data_name '_variogram.1D']);
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

c = corr(y,x,'Type',CorrType,'Rows', 'pairwise')
cn = corr(Data1_null,x,'Type',CorrType,'Rows', 'pairwise');
p = sum(abs(c)<abs(cn))/20000

if ifDeOutlier
    mask = subfun_FeatureSpace_DeOutlier(x);
    y = y(mask>0);
    x = x(mask>0);
end

figure('Color',[1 1 1],'Position',[0 0 300 300],'Units','pixels');
hold on;
% subfun_plot_scatter_density_isopotential( x , y , valueorrank , ifdensity , ifisopotential , ifregressline )
subfun_plot_scatter_regression(y, x, [-12 9] , [2.5 4.5]);



%% compare CBF-score map and previous scaling-PNC map (Fig S2e)
x = all_other_maps.maps{60};
Data_name = 'CBF_PC1_Val_L12_bil';  % 'CBF_PC1_Val_bil', 'CBF_Mean_PScs_bil'
y = load([Group_dir Data_name '.1D']);
Data1_null = load([Group_dir Data_name '_variogram.1D']);
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

c = corr(y,x,'Type',CorrType,'Rows', 'pairwise')
cn = corr(Data1_null,x,'Type',CorrType,'Rows', 'pairwise');
p = sum(abs(c)<abs(cn))/20000

if ifDeOutlier
    mask = subfun_FeatureSpace_DeOutlier(x);
    y = y(mask>0);
    x = x(mask>0);
end

figure('Color',[1 1 1],'Position',[0 0 300 300],'Units','pixels');
hold on;
% subfun_plot_scatter_density_isopotential( x , y , valueorrank , ifdensity , ifisopotential , ifregressline )
subfun_plot_scatter_regression(y, x, [-12 9] , [0.6 1.4]);



%% write to SourceData EDFig2
Data_name = 'CBF_PC1_Val_L12_bil';  % 'CBF_PC1_Val_bil', 'CBF_Mean_PScs_bil'
y = load([Group_dir Data_name '.1D']);
outSD = [y all_other_maps.maps{54} all_other_maps.maps{64} all_other_maps.maps{56} all_other_maps.maps{57} all_other_maps.maps{75} all_other_maps.maps{60}];
for i = 2:7
    x = outSD(:,i);
    mask = CBFsnr_mean_bil>0 & abs(x-mean(x))<3*std(x);
    outSD(~mask,i) = nan;
end




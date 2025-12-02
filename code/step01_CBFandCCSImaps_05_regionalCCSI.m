clc
clear all
close all
set(0,'defaultfigurecolor',[1 1 1]);% set figure background color
addpath(genpath('/Users/guofanhua/Desktop/gfh/tools/MatlabTools/BrainSpace-latest/matlab'))

%% ======================= compare regional CCSI

%% input data
step01_CBFandCCSImaps_00_init;
step01_CBFandCCSImaps_00_ComputeBasicResults
CorrType = 'Pearson';  % Pearson Spearman


%% regional CCSI without SNR threshold (Fig S3a)
[~,p] = ttest(Data_PScs_uni');
BHp = fdr_BH(p,0.05)';
SNRmean = mean(subfun_merge_bil2uni(Data_SNR')',2);
SNRstd = std(subfun_merge_bil2uni(Data_SNR'))';

tmp = atlas_info.Lobe(1:180);
atlas_lobeindex = zeros(size(tmp));
for i = 1:numel(atlas_lobeindex)
    switch tmp{i}
        case 'Fr'
            atlas_lobeindex(i) = 1;
        case 'Ins'
            atlas_lobeindex(i) = 2;
        case 'Occ'
            atlas_lobeindex(i) = 3;
        case 'Par'
            atlas_lobeindex(i) = 4;
        case 'Temp'
            atlas_lobeindex(i) = 5;
    end
end

tmp = table(mean(Data_PScs_uni,2), std(Data_PScs_uni')'/sqrt(size(Data_PScs_uni,2)), atlas_info.Lobe(1:180), atlas_lobeindex, p',...
            BHp, SNRmean, SNRstd, 'VariableNames',...
            {'mean', 'ste', 'lobe', 'lobeIndex', 'p',...
            'BHp', 'SNRmean', 'SNRstd'});
tmp = sortrows(tmp,6,'ascend');
tmp = sortrows(tmp,3);


figure('Color',[1 1 1],'Position',[0 0 1200 400],'Units','pixels');
hold on;
b = bar(sum(tmp.lobeIndex<1)+1:sum(tmp.lobeIndex<2), tmp.mean(sum(tmp.lobeIndex<1)+1:sum(tmp.lobeIndex<2)), 'FaceColor', 'white', 'EdgeColor', [0, 200, 180]/255, 'LineWidth', 1);
b = bar(sum(tmp.lobeIndex<2)+1:sum(tmp.lobeIndex<3), tmp.mean(sum(tmp.lobeIndex<2)+1:sum(tmp.lobeIndex<3)), 'FaceColor', 'white', 'EdgeColor', [242, 142, 43]/255, 'LineWidth', 1);
b = bar(sum(tmp.lobeIndex<3)+1:sum(tmp.lobeIndex<4), tmp.mean(sum(tmp.lobeIndex<3)+1:sum(tmp.lobeIndex<4)), 'FaceColor', 'white', 'EdgeColor', [89, 161, 79]/255, 'LineWidth', 1);
b = bar(sum(tmp.lobeIndex<4)+1:sum(tmp.lobeIndex<5), tmp.mean(sum(tmp.lobeIndex<4)+1:sum(tmp.lobeIndex<5)), 'FaceColor', 'white', 'EdgeColor', [225, 87, 89]/255, 'LineWidth', 1);
b = bar(sum(tmp.lobeIndex<5)+1:sum(tmp.lobeIndex<6), tmp.mean(sum(tmp.lobeIndex<5)+1:sum(tmp.lobeIndex<6)), 'FaceColor', 'white', 'EdgeColor', [150, 103, 185]/255, 'LineWidth', 1);
errorbar(1:180,tmp.mean,tmp.ste,'k.','LineWidth', 0.5);
legend({'FL','Ins','OL','PL','TL'})
set(gca, 'XColor', 'none', 'xTick', []);
set(gca, 'ylim', [-pi/4 pi/2], 'yTick', [-pi/4 0 pi/2], 'ytickLabel',{'-π/4','0','π/2'});


% out stats
tmplabel = {'FL','Ins','OL','PL','TL'};
for i = 1:numel(tmplabel)
    fprintf([tmplabel{i} ' spin-p-FDR < 0.05 have ' num2str(sum(tmp.BHp(sum(tmp.lobeIndex<i)+1:sum(tmp.lobeIndex<i+1))<0.05)) '/' num2str(sum(tmp.lobeIndex==i)) '\n'])
end



%% find the threshold of SNR-CCSI
x = mean(Data_SNR,2);
Data_name = 'CBF_Mean_PScs_bil';  % 'CBF_PC1_Val_bil', 'CBF_Mean_PScs_bil'
y = load([Group_dir Data_name '.1D']);
Data1_null = load([Group_dir Data_name '_null.1D']);


% Spearman
C_trend_up = [];
p_trend_up = [];
C_trend_low = [];
p_trend_low = [];
for i = 2:0.1:7
    mask = x>=i;
    ty = y(mask);
    tx = x(mask);
    tnull = Data1_null(mask,:);
    c = corr(ty,tx,'Type',CorrType);
    cn = corr(tnull,tx,'Type',CorrType);
    p = sum(abs(c)<abs(cn))/20000;
    C_trend_up = [C_trend_up; c];
    p_trend_up = [p_trend_up; p];
    
    mask = x<i;
    ty = y(mask);
    tx = x(mask);
    tnull = Data1_null(mask,:);
    c = corr(ty,tx,'Type',CorrType);
    cn = corr(tnull,tx,'Type',CorrType);
    p = sum(abs(c)<abs(cn))/20000;
    C_trend_low = [C_trend_low; c];
    p_trend_low = [p_trend_low; p];
end




%% compare CBF-score map and previous PET-CBF map (Fig S3cde)
% whole SNR 3c
x = CBFsnr_mean_bil;
Data_name = 'CBF_Mean_PScs_bil';  % 'CBF_PC1_Val_bil', 'CBF_Mean_PScs_bil'
y = load([Group_dir Data_name '.1D']);
Data1_null = load([Group_dir Data_name '_null.1D']);
ifdensity = 1;
ifisopotential = 1;
ifregressline = 0;
valueorrank = 2;  % 1:value 2:rank
ifDeOutlier = 0;

mask = CBFsnr_mean_bil>0;
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


% SNR<4 3d
x = CBFsnr_mean_bil;
Data_name = 'CBF_Mean_PScs_bil';  % 'CBF_PC1_Val_bil', 'CBF_Mean_PScs_bil'
y = load([Group_dir Data_name '.1D']);
Data1_null = load([Group_dir Data_name '_null.1D']);
ifdensity = 1;
ifisopotential = 1;
ifregressline = 0;
valueorrank = 2;  % 1:value 2:rank
ifDeOutlier = 0;

mask = CBFsnr_mean_bil<4;
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


% SNR>4 3e
x = CBFsnr_mean_bil;
Data_name = 'CBF_Mean_PScs_bil';  % 'CBF_PC1_Val_bil', 'CBF_Mean_PScs_bil'
y = load([Group_dir Data_name '.1D']);
Data1_null = load([Group_dir Data_name '_null.1D']);
ifdensity = 1;
ifisopotential = 1;
ifregressline = 0;
valueorrank = 2;  % 1:value 2:rank
ifDeOutlier = 0;

mask = CBFsnr_mean_bil>=4;
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





%% regional CCSI with SNR>4 threshold (Fig 1e)
[~,p] = ttest(Data_PScs_uni');
BHp = fdr_BH(p,0.05)';
SNRmean = mean(subfun_merge_bil2uni(Data_SNR')',2);
SNRstd = std(subfun_merge_bil2uni(Data_SNR'))';

tmp = atlas_info.Lobe(1:180);
atlas_lobeindex = zeros(size(tmp));
for i = 1:numel(atlas_lobeindex)
    switch tmp{i}
        case 'Fr'
            atlas_lobeindex(i) = 1;
        case 'Ins'
            atlas_lobeindex(i) = 2;
        case 'Occ'
            atlas_lobeindex(i) = 3;
        case 'Par'
            atlas_lobeindex(i) = 4;
        case 'Temp'
            atlas_lobeindex(i) = 5;
    end
end

tmp = table(mean(Data_PScs_uni,2), std(Data_PScs_uni')'/sqrt(size(Data_PScs_uni,2)), atlas_info.Lobe(1:180), atlas_lobeindex, p',...
            BHp, SNRmean, SNRstd, atlas_info.regionName(1:180), 'VariableNames',...
            {'mean', 'ste', 'lobe', 'lobeIndex', 'p',...
            'BHp', 'SNRmean', 'SNRstd', 'regionName'});
tmp = sortrows(tmp,6,'ascend');
tmp = sortrows(tmp,3);
tmp(tmp.SNRmean<4,:) = [];

figure('Color',[1 1 1],'Position',[0 0 1200 400],'Units','pixels');
hold on;
b = bar(sum(tmp.lobeIndex<1)+1:sum(tmp.lobeIndex<2), tmp.mean(sum(tmp.lobeIndex<1)+1:sum(tmp.lobeIndex<2)), 'FaceColor', 'white', 'EdgeColor', [0, 200, 180]/255, 'LineWidth', 1);
b = bar(sum(tmp.lobeIndex<2)+1:sum(tmp.lobeIndex<3), tmp.mean(sum(tmp.lobeIndex<2)+1:sum(tmp.lobeIndex<3)), 'FaceColor', 'white', 'EdgeColor', [242, 142, 43]/255, 'LineWidth', 1);
b = bar(sum(tmp.lobeIndex<3)+1:sum(tmp.lobeIndex<4), tmp.mean(sum(tmp.lobeIndex<3)+1:sum(tmp.lobeIndex<4)), 'FaceColor', 'white', 'EdgeColor', [89, 161, 79]/255, 'LineWidth', 1);
b = bar(sum(tmp.lobeIndex<4)+1:sum(tmp.lobeIndex<5), tmp.mean(sum(tmp.lobeIndex<4)+1:sum(tmp.lobeIndex<5)), 'FaceColor', 'white', 'EdgeColor', [225, 87, 89]/255, 'LineWidth', 1);
b = bar(sum(tmp.lobeIndex<5)+1:sum(tmp.lobeIndex<6), tmp.mean(sum(tmp.lobeIndex<5)+1:sum(tmp.lobeIndex<6)), 'FaceColor', 'white', 'EdgeColor', [150, 103, 185]/255, 'LineWidth', 1);
errorbar(1:numel(tmp.mean),tmp.mean,tmp.ste,'k.','LineWidth', 0.5);
legend({'FL','Ins','OL','PL','TL'})
set(gca, 'XColor', 'none', 'xTick', []);
set(gca, 'ylim', [-pi/4 pi/2], 'yTick', [-pi/4 0 pi/2], 'ytickLabel',{'-π/4','0','π/2'});


% out stats
tmplabel = {'FL','Ins','OL','PL','TL'};
for i = 1:numel(tmplabel)
    fprintf([tmplabel{i} ' spin-p-FDR < 0.05 have ' num2str(sum(tmp.BHp(sum(tmp.lobeIndex<i)+1:sum(tmp.lobeIndex<i+1))<0.05)) '/' num2str(sum(tmp.lobeIndex==i)) '\n'])
end


%% one-way repeated measures ANOVA
clc
data = Data_PScs_uni(CBFsnr_mean_uni>4,:)';
%----- 输入：data 为 [Nsub × Nroi] -----
[Nsub, Nroi] = size(data);
roiNames = compose("ROI_%03d", 1:Nroi);   % 生成变量名 ROI_001 ... ROI_150
T = array2table(data, "VariableNames", cellstr(roiNames));
T.Subject = (1:Nsub)';                     % 被试编号

%----- 构造重复测量设计（“Condition”即脑区）-----
within = table(categorical(roiNames'),'VariableNames',{'Condition'});

%----- 建模：A-Z 这种区间写法需用首尾变量名拼出公式 -----
formula = sprintf('%s-%s ~ 1', T.Properties.VariableNames{1}, ...
                              T.Properties.VariableNames{Nroi});
rm = fitrm(T, formula, "WithinDesign", within);

%----- 单因素重复测量 ANOVA（含球形性校正列）-----
ranovatbl = ranova(rm, "WithinModel", "Condition");
disp(ranovatbl)

% 可选：球形性检验与校正因子
mauchlyTbl = mauchly(rm);   % Mauchly’s test
epsTbl     = epsilon(rm);   % GG/HF epsilon
disp(mauchlyTbl), disp(epsTbl)


% 找到Condition行和对应的误差行
ss_effect = ranovatbl.SumSq(3);
ss_error  = ranovatbl.SumSq(4);

% 计算偏η²
eta_p2 = ss_effect / (ss_effect + ss_error);

fprintf("Partial eta squared (η_p^2) for Condition = %.4f\n", eta_p2);




%% write to map plot by R
OutNames = {'CBF_Mean_PScs_qc_bil.csv'};
SNRmean = CBFsnr_mean_bil;
OutDatas = {CBF_meanPScs_bil(SNRmean>4)};
Datalabel = [1:360]';
Datalabel(SNRmean<=4) = [];

for i = 1:numel(OutNames)
    tmp = subfun_write_to_MapPlot_byR(OutDatas{i},Datalabel,atlas_info,[Group_dir OutNames{i}]);
end






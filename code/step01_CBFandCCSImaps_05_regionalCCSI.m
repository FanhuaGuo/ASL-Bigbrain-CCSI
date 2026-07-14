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
step01_CBFandCCSImaps_00_ComputeBasicResults
CorrType = 'Pearson';  % Pearson Spearman


%% ======================= compare left and right hemi (symmetric analysis)
%% compare average lh and rh  (Fig S3bc)
CCSI_lh = Data_PScs_bil(1:180,:);
CCSI_rh = Data_PScs_bil(181:360,:);
CCSI_lh_mean = mean(CCSI_lh(CBFsnr_mean_uni>4,:))';
CCSI_rh_mean = mean(CCSI_rh(CBFsnr_mean_uni>4,:))';
CCSI_lh_std = std(CCSI_lh(CBFsnr_mean_uni>4,:))';
CCSI_rh_std = std(CCSI_rh(CBFsnr_mean_uni>4,:))';

% test 
[~,p_mean] = ttest(CCSI_lh_mean,CCSI_rh_mean)
[~,p_std] = ttest(CCSI_lh_std,CCSI_rh_std)


% plot
figure('Color',[1 1 1],'Position',[0 0 300 300],'Units','pixels');
subplot(1,2,1);
hold on;
% 画每个被试的连线
for i = 1:numel(CCSI_lh_mean)
    plot([1 2], [CCSI_lh_mean(i) CCSI_rh_mean(i)], '-o', ...
        'Color', [0.7 0.7 0.7], 'MarkerSize', 5, 'Linewidth', 0.5);
end
% 画均值（加粗）
plot([1 2], [mean(CCSI_lh_mean) mean(CCSI_rh_mean)], '-o', ...
    'Color', 'k', 'LineWidth', 2, 'MarkerSize', 8, 'Linewidth', 1.5);
xlim([0.5 2.5])
xticks([1 2])
xticklabels({'Left','Right'})
ylim([0 1]);
yticks([0 1]);
ylabel('CCSI')
title('mean')
box off


subplot(1,2,2);
hold on;
% 画每个被试的连线
for i = 1:numel(CCSI_lh_std)
    plot([1 2], [CCSI_lh_std(i) CCSI_rh_std(i)], '-o', ...
        'Color', [0.7 0.7 0.7], 'MarkerSize', 5, 'Linewidth', 0.5);
end
% 画均值（加粗）
plot([1 2], [mean(CCSI_lh_std) mean(CCSI_rh_std)], '-o', ...
    'Color', 'k', 'LineWidth', 2, 'MarkerSize', 8, 'Linewidth', 1.5);
xlim([0.5 2.5])
xticks([1 2])
xticklabels({'Left','Right'})
ylim([0 1]);
yticks([0 1]);
ylabel('CCSI')
title('std')
box off


% write to SourceData EDFig3b
outSD = [CCSI_lh_mean CCSI_lh_std CCSI_rh_mean CCSI_rh_std];


% plot histogram
data = mean(CCSI_lh,2);
[h,p] = kstest((data-mean(data))/std(data))
data = mean(CCSI_rh,2);
[h,p] = kstest((data-mean(data))/std(data))

data = [mean(CCSI_lh,2); mean(CCSI_rh,2)];
group = [ones(size(mean(CCSI_lh,2))); 2*ones(size(mean(CCSI_rh,2)))];
p_var = vartestn(data, group, 'TestType','Bartlett');



%% spatial distribution compare (asymmetric)
Data_name = ['CBF_Mean_PScs_uncorr_L12_bil'];
y = load([Group_dir Data_name '.1D']);
Data1_null = load([Group_dir Data_name '_variogram.1D']);
% CorrType = 'Pearson';  % Pearson Spearman

% plot
figure('Color',[1 1 1],'Position',[0 0 300 300],'Units','pixels');
hold on;
subfun_plot_scatter_regression(y(1:end/2) , y(end/2+1:end) , [-0.5 1.1] , [-0.5 1.1]);

% variogram p
c = corr(y(1:end/2), y(end/2+1:end),'Type',CorrType,'Rows', 'pairwise')
cn = corr(Data1_null(1:end/2,:), y(end/2+1:end),'Type',CorrType,'Rows', 'pairwise');
p = sum(abs(c)<abs(cn))/20000


% write to SourceData EDFig3b
outSD = [y(1:end/2) y(end/2+1:end)];



%% ======================= regional CCSI and QC
%% single region, single subject --- compare laminar profiles of CBF/qT1/CSI (Fig 3b)
region_label = 1;
subject_label = 1;

tmpCBF = Data_scan_bil{subject_label}(region_label,:);
tmpqT1 = Data_T1v_bil{subject_label}(region_label,:);
tmpCSI = BB_LaminarProfile(region_label,:);

tmpCBF = zscore(tmpCBF);
tmpqT1 = zscore(tmpqT1);
tmpCSI = zscore(tmpCSI);

% plot
figure('Color',[1 1 1],'Position',[0 0 400 300],'Units','pixels');
hold on;
plot(tmpqT1,'b-','linewidth',1);
plot(tmpCBF,'r-','linewidth',1);
plot(tmpCSI,'k-','linewidth',1);

xlim([0 numel(tmpCBF)+1])
xticks([])
ylim([-3 3]);
yticks([]);
box on


% write to SourceData EDFig3b
outSD = [tmpCSI' tmpCBF' tmpqT1'];



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
b = bar(sum(tmp.lobeIndex<1)+1:sum(tmp.lobeIndex<2), tmp.mean(sum(tmp.lobeIndex<1)+1:sum(tmp.lobeIndex<2)), 'FaceColor', 'white', 'EdgeColor', [0, 200, 180]/255, 'LineWidth', 0.5);
b = bar(sum(tmp.lobeIndex<2)+1:sum(tmp.lobeIndex<3), tmp.mean(sum(tmp.lobeIndex<2)+1:sum(tmp.lobeIndex<3)), 'FaceColor', 'white', 'EdgeColor', [242, 142, 43]/255, 'LineWidth', 0.5);
b = bar(sum(tmp.lobeIndex<3)+1:sum(tmp.lobeIndex<4), tmp.mean(sum(tmp.lobeIndex<3)+1:sum(tmp.lobeIndex<4)), 'FaceColor', 'white', 'EdgeColor', [89, 161, 79]/255, 'LineWidth', 0.5);
b = bar(sum(tmp.lobeIndex<4)+1:sum(tmp.lobeIndex<5), tmp.mean(sum(tmp.lobeIndex<4)+1:sum(tmp.lobeIndex<5)), 'FaceColor', 'white', 'EdgeColor', [225, 87, 89]/255, 'LineWidth', 0.5);
b = bar(sum(tmp.lobeIndex<5)+1:sum(tmp.lobeIndex<6), tmp.mean(sum(tmp.lobeIndex<5)+1:sum(tmp.lobeIndex<6)), 'FaceColor', 'white', 'EdgeColor', [150, 103, 185]/255, 'LineWidth', 0.5);
errorbar(1:numel(tmp.mean),tmp.mean,tmp.ste,'k.','LineWidth', 0.5);
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
Data_name = ['CBF_Mean_PScs_L' num2str(num_layer) '_bil'];  % 'CBF_PC1_Val_bil', 'CBF_Mean_PScs_bil'
NullType = 'variogram';    % variogram  spin  eigenstrapping
y = load([Group_dir Data_name '.1D']);
Data1_null = load([Group_dir Data_name '_' NullType '.1D']);


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
    c = corr(ty,tx,'Type',CorrType, 'Rows', 'pairwise');
    cn = corr(tnull,tx,'Type',CorrType, 'Rows', 'pairwise');
    p = sum(abs(c)<abs(cn))/20000;
    C_trend_up = [C_trend_up; c];
    p_trend_up = [p_trend_up; p];
    
    mask = x<i;
    ty = y(mask);
    tx = x(mask);
    tnull = Data1_null(mask,:);
    c = corr(ty,tx,'Type',CorrType, 'Rows', 'pairwise');
    cn = corr(tnull,tx,'Type',CorrType, 'Rows', 'pairwise');
    p = sum(abs(c)<abs(cn))/20000;
    C_trend_low = [C_trend_low; c];
    p_trend_low = [p_trend_low; p];
end




%% compare CBF-score map and previous PET-CBF map (Fig S3cde)
NullType = 'variogram';    % variogram  spin  eigenstrapping

% whole SNR 3g
x = CBFsnr_mean_bil;
Data_name = ['CBF_Mean_PScs_L' num2str(num_layer) '_bil'];  % 'CBF_PC1_Val_bil', 'CBF_Mean_PScs_bil'
y = load([Group_dir Data_name '.1D']);
Data1_null = load([Group_dir Data_name '_' NullType '.1D']);
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
% subfun_plot_scatter_density_isopotential( x , y , valueorrank , ifdensity , ifisopotential , ifregressline )
subfun_plot_scatter_regression(x , y , [0 10] , [-0.5 1.2]);


% SNR<4 3h
x = CBFsnr_mean_bil;
Data_name = ['CBF_Mean_PScs_L' num2str(num_layer) '_bil'];  % 'CBF_PC1_Val_bil', 'CBF_Mean_PScs_bil'
y = load([Group_dir Data_name '.1D']);
Data1_null = load([Group_dir Data_name '_' NullType '.1D']);
ifdensity = 1;
ifisopotential = 1;
ifregressline = 0;
valueorrank = 2;  % 1:value 2:rank
ifDeOutlier = 0;

mask = CBFsnr_mean_bil<4;
x = x(mask);
y = y(mask);
Data1_null = Data1_null(mask,:);

c = corr(y,x,'Type',CorrType, 'Rows', 'pairwise')
cn = corr(Data1_null,x,'Type',CorrType, 'Rows', 'pairwise');
p = sum(abs(c)<abs(cn))/20000

if ifDeOutlier
    mask = subfun_FeatureSpace_DeOutlier(x);
    y = y(mask>0);
    x = x(mask>0);
end

figure('Color',[1 1 1],'Position',[0 0 300 300],'Units','pixels');
hold on;
% subfun_plot_scatter_density_isopotential( x , y , valueorrank , ifdensity , ifisopotential , ifregressline )
subfun_plot_scatter_regression(x , y , [0 10] , [-0.5 1.2]);


% SNR>4 3i
x = CBFsnr_mean_bil;
Data_name = ['CBF_Mean_PScs_L' num2str(num_layer) '_bil'];  % 'CBF_PC1_Val_bil', 'CBF_Mean_PScs_bil'
y = load([Group_dir Data_name '.1D']);
Data1_null = load([Group_dir Data_name '_' NullType '.1D']);
ifdensity = 1;
ifisopotential = 1;
ifregressline = 0;
valueorrank = 2;  % 1:value 2:rank
ifDeOutlier = 0;

mask = CBFsnr_mean_bil>=4;
x = x(mask);
y = y(mask);
Data1_null = Data1_null(mask,:);

c = corr(y,x,'Type',CorrType, 'Rows', 'pairwise')
cn = corr(Data1_null,x,'Type',CorrType, 'Rows', 'pairwise');
p = sum(abs(c)<abs(cn))/20000

if ifDeOutlier
    mask = subfun_FeatureSpace_DeOutlier(x);
    y = y(mask>0);
    x = x(mask>0);
end

figure('Color',[1 1 1],'Position',[0 0 300 300],'Units','pixels');
hold on;
% subfun_plot_scatter_density_isopotential( x , y , valueorrank , ifdensity , ifisopotential , ifregressline )
subfun_plot_scatter_regression(x , y , [0 10] , [-0.5 1.2]);


% write to SourceData EDFig3ghi
y = load([Group_dir Data_name '.1D']);
outSD = [CBFsnr_mean_bil y];


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
b = bar(sum(tmp.lobeIndex<1)+1:sum(tmp.lobeIndex<2), tmp.mean(sum(tmp.lobeIndex<1)+1:sum(tmp.lobeIndex<2)), 'FaceColor', 'white', 'EdgeColor', [0, 200, 180]/255, 'LineWidth', 0.5);
b = bar(sum(tmp.lobeIndex<2)+1:sum(tmp.lobeIndex<3), tmp.mean(sum(tmp.lobeIndex<2)+1:sum(tmp.lobeIndex<3)), 'FaceColor', 'white', 'EdgeColor', [242, 142, 43]/255, 'LineWidth', 0.5);
b = bar(sum(tmp.lobeIndex<3)+1:sum(tmp.lobeIndex<4), tmp.mean(sum(tmp.lobeIndex<3)+1:sum(tmp.lobeIndex<4)), 'FaceColor', 'white', 'EdgeColor', [89, 161, 79]/255, 'LineWidth', 0.5);
b = bar(sum(tmp.lobeIndex<4)+1:sum(tmp.lobeIndex<5), tmp.mean(sum(tmp.lobeIndex<4)+1:sum(tmp.lobeIndex<5)), 'FaceColor', 'white', 'EdgeColor', [225, 87, 89]/255, 'LineWidth', 0.5);
b = bar(sum(tmp.lobeIndex<5)+1:sum(tmp.lobeIndex<6), tmp.mean(sum(tmp.lobeIndex<5)+1:sum(tmp.lobeIndex<6)), 'FaceColor', 'white', 'EdgeColor', [150, 103, 185]/255, 'LineWidth', 0.5);
errorbar(1:numel(tmp.mean),tmp.mean,tmp.ste,'k.','LineWidth', 0.5);
legend({'FL','Ins','OL','PL','TL'})
set(gca, 'XColor', 'none', 'xTick', []);
set(gca, 'ylim', [-pi/4 pi/2], 'yTick', [-pi/4 0 pi/2], 'ytickLabel',{'-π/4','0','π/2'});


% out stats
tmplabel = {'FL','Ins','OL','PL','TL'};
for i = 1:numel(tmplabel)
    fprintf([tmplabel{i} ' spin-p-FDR < 0.05 have ' num2str(sum(tmp.BHp(sum(tmp.lobeIndex<i)+1:sum(tmp.lobeIndex<i+1))<0.05)) '/' num2str(sum(tmp.lobeIndex==i)) '\n'])
end


% write to SourceData Fig3d
outSD = table(atlas_info.region(1:180),Data_PScs_uni);
outSD(SNRmean<4,:) = [];

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



%% show laminar profile of CBF/CSI/qT1 in different regions
close all
% find label
[~,p] = ttest(Data_PScs_uni');
BHp = fdr_BH(p,0.05)';
tmp = [mean(Data_PScs_uni,2), [1:180]', BHp];
tmp = tmp(CBFsnr_mean_uni>4,:);
CCSI_rank = sortrows(tmp,1,'descend');
% pick_regions = [4 1 53 153 133];
pick_regions = [4 1 50 139 107 167];
pick_rank = [];
for i = 1:numel(pick_regions)
    tmp = find(CCSI_rank(:,2)==pick_regions(i));
    pick_rank = [pick_rank; tmp];
end
pick_rank = 100*pick_rank/size(CCSI_rank,1);


figure('Color',[1 1 1],'Position',[0 0 600 600],'Units','pixels');
for ima = 1:6
    tt = pick_regions(ima);
    subplot(2,3,ima); hold on;
    
    tmp1 = [];
    tmp2 = [];
    for i = 1:numel(Data_scan_uni)
        tmp1 = [tmp1; zscore(Data_scan_uni{i}(tt,1:end-1))];
    %     plot(zscore(Data_scan_uni{i}(tt,1:end-1)),'r-','linewidth',0.5);
        tmp2 = [tmp2; zscore(Data_T1v_uni{i}(tt,1:end-1))];
    end
    tmpCBF = mean(tmp1,1);
    tmpCBFste = std(tmp1);%/sqrt(numel(Data_scan_uni));
    tmpqT1 = mean(tmp2,1);
    tmpqT1ste = std(tmp2);%/sqrt(numel(Data_scan_uni));
    tmpCSI = zscore(Data_BBLP_uni(tt,1:end-1))';
    
    fill([[1:num_layer-1]  flip([1:num_layer-1])]', [tmpqT1+tmpqT1ste flip(tmpqT1-tmpqT1ste)]', ...
         [0.7 0.7 1], 'EdgeColor', 'none');    
    fill([[1:num_layer-1]  flip([1:num_layer-1])]', [tmpCBF+tmpCBFste flip(tmpCBF-tmpCBFste)]', ...
         [1 0.7 0.7], 'EdgeColor', 'none');    
    plot(tmpqT1,'b-','linewidth',1);
    plot(tmpCBF,'r-','linewidth',1);
    plot(tmpCSI,'k-','linewidth',1);

    xlim([0 numel(tmpCBF)+2])
    xticks([])
    ylim([-3 3]);
    yticks([]);
    box on
end


% write to SourceData Fig3f
outSD = {};
for ima = 1:6
    tt = pick_regions(ima);
    
    tmp1 = [];
    tmp2 = [];
    for i = 1:numel(Data_scan_uni)
        tmp1 = [tmp1; zscore(Data_scan_uni{i}(tt,1:end-1))];
    %     plot(zscore(Data_scan_uni{i}(tt,1:end-1)),'r-','linewidth',0.5);
        tmp2 = [tmp2; zscore(Data_T1v_uni{i}(tt,1:end-1))];
    end
    tmpCBF = mean(tmp1,1);
    tmpCBFste = std(tmp1);%/sqrt(numel(Data_scan_uni));
    tmpqT1 = mean(tmp2,1);
    tmpqT1ste = std(tmp2);%/sqrt(numel(Data_scan_uni));
    tmpCSI = zscore(Data_BBLP_uni(tt,1:end-1))';
    
    outSD{ima} = [tmpCSI tmpCBF' tmpCBFste' tmp1' tmpqT1' tmpqT1ste' tmp2'];
end

%% write to map plot by R
OutNames = {['CBF_Mean_PScs_qc_bil_L' num2str(num_layer) '.csv'],...
            ['T1vBB_Mean_PScs_qc_bil_L' num2str(num_layer) '.csv']};
SNRmean = CBFsnr_mean_bil;
OutDatas = {CBF_meanPScs_bil_correct(SNRmean>4),...
            T1vBB_meanPScs_bil_correct(SNRmean>4)};
Datalabel = [1:360]';
Datalabel(SNRmean<=4) = [];

for i = 1:numel(OutNames)
    tmp = subfun_write_to_MapPlot_byR(OutDatas{i},Datalabel,atlas_info,[Group_dir OutNames{i}]);
end






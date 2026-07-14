clc
clear all
close all
set(0,'defaultfigurecolor',[1 1 1]);% set figure background color
addpath('/Users/guofanhua/Desktop/gfh/tools/MatlabTools/BrainSpace-latest/matlab')

%% ======================= step 1 for disposal all data

%% input data
step05_retest_00_init;
step01_CBFandCCSImaps_00_init;

%% compute cosine similarity between CBF and Cell-body staining intensity (CCSI)
% s1
Data1_PScs_bil = subfun_ProfileSimilarity(Data1_scan_bil,Data_BBLP_bil);
Data1_PScs_uni = subfun_ProfileSimilarity(Data1_scan_uni,Data_BBLP_uni);
Data1_PScs_bil = pi-acos(Data1_PScs_bil) -pi/2;
Data1_PScs_uni = pi-acos(Data1_PScs_uni) -pi/2;

% s2
Data2_PScs_bil = subfun_ProfileSimilarity(Data2_scan_bil,Data_BBLP_bil);
Data2_PScs_uni = subfun_ProfileSimilarity(Data2_scan_uni,Data_BBLP_uni);
Data2_PScs_bil = pi-acos(Data2_PScs_bil) -pi/2;
Data2_PScs_uni = pi-acos(Data2_PScs_uni) -pi/2;


%% mean or PCA combine all subjects into 1 map
tmp = subfun_PCAtoPC(zscore(Data1_Val_uni));
CBF1_ValPC1_uni = tmp(:,1);
[tmp, explained] = subfun_PCAtoPC(zscore(Data1_Val_bil));
CBF1_ValPC1_bil = tmp(:,1);
tmp = subfun_PCAtoPC(zscore(Data2_Val_uni));
CBF2_ValPC1_uni = tmp(:,1);
[tmp, explained] = subfun_PCAtoPC(zscore(Data2_Val_bil));
CBF2_ValPC1_bil = tmp(:,1);


CBF1_meanPScs_bil = mean(Data1_PScs_bil,2);
CBF1_meanPScs_uni = mean(Data1_PScs_uni,2);
CBF2_meanPScs_bil = mean(Data2_PScs_bil,2);
CBF2_meanPScs_uni = mean(Data2_PScs_uni,2);

CBFsnr_mean_bil = mean(Data_SNR,2);
CBFsnr_mean_uni = subfun_merge_bil2uni(CBFsnr_mean_bil);



%%
icc_ccsi = [];
icc_cbf = [];
for i = 1:360
    icc_ccsi = [icc_ccsi; subfun_compute_icc3_1([Data1_PScs_bil(i,:)' Data2_PScs_bil(i,:)'])];
    icc_cbf = [icc_cbf; subfun_compute_icc3_1([Data1_Val_bil(i,:)' Data2_Val_bil(i,:)'])];
end
% icc = icc(CBFsnr_mean_bil>4);



%% Bland–Altman
% 假设你的数据是两个向量:
test = CBF1_meanPScs_bil(CBFsnr_mean_bil>4);     %[nROI x 1]  (test session 的 CCSI)
retest = CBF2_meanPScs_bil(CBFsnr_mean_bil>4);   %[nROI x 1] (retest session 的 CCSI)

% 1. 计算均值和差值
meanVals = (test + retest) / 2;   % 横轴
diffVals = test - retest;         % 纵轴

% 2. 计算 bias 和一致性限
bias = mean(diffVals);            
sd_diff = std(diffVals);          
loa_upper = bias + 1.96 * sd_diff; % 上一致性限
loa_lower = bias - 1.96 * sd_diff; % 下一致性限

% 3. 作图
figure; hold on;
scatter(meanVals, diffVals, 40, 'k', 'filled'); % 散点图
yline(bias, 'r-', 'LineWidth', 2, 'Label', 'Bias'); 
yline(loa_upper, 'b--', 'LineWidth', 2, 'Label', '+1.96 SD');
yline(loa_lower, 'b--', 'LineWidth', 2, 'Label', '-1.96 SD');
xlabel('Mean of Test and Retest CCSI');
ylabel('Difference (Test - Retest)');
title('Bland-Altman Plot for Test–Retest CCSI');
grid on;
set(gca, 'FontSize', 12);

% 4. 输出结果
fprintf('Bias = %.4f\n', bias);
fprintf('Upper LoA = %.4f\n', loa_upper);
fprintf('Lower LoA = %.4f\n', loa_lower);



%% roi-wise icc
% 假设你的数据矩阵如下：
% test_data  = [nSubjects x nROIs]   % test session 的 CCSI
% retest_data = [nSubjects x nROIs]  % retest session 的 CCSI
% nSubjects = 14, nROIs = 360 (HCP-MMP1 parcellation)
test_data = Data1_PScs_bil(CBFsnr_mean_bil>4,:)';     %[nROI x 1]  (test session 的 CCSI)
retest_data = Data2_PScs_bil(CBFsnr_mean_bil>4,:)';   %[nROI x 1] (retest session 的 CCSI)

nSubjects = size(test_data,1);
nROIs = size(test_data,2);

icc_vals = nan(nROIs,1);

% 逐 ROI 计算 ICC(2,1) - 双向随机效应，单次测量，绝对一致性
for roi = 1:nROIs
    Y = [test_data(:,roi), retest_data(:,roi)];  % [nSubjects x 2]
    
    % 去掉 NaN
    Y(any(isnan(Y),2),:) = [];
    
    % 计算方差分解
    n = size(Y,1); % 有效样本数
    k = size(Y,2); % 重复次数 (这里是2: test & retest)
    mean_per_subject = mean(Y,2);
    mean_per_rater   = mean(Y,1);
    grand_mean = mean(Y(:));
    
    % 平方和
    SS_between = k * sum((mean_per_subject - grand_mean).^2);
    SS_within  = sum(sum((Y - mean_per_subject).^2));
    SS_rater   = n * sum((mean_per_rater - grand_mean).^2);
    
    % 均方
    MS_between = SS_between / (n-1);
    MS_within  = SS_within / ((n-1)*(k-1));
    MS_rater   = SS_rater / (k-1);
    
    % ICC(2,1) 公式
    icc = (MS_between - MS_within) / ...
          (MS_between + (k-1)*MS_within + (k/n)*(MS_rater - MS_within));
    
    icc_vals(roi) = icc;
end

% 输出结果
fprintf('ROI-wise ICC summary:\n');
fprintf('Mean ICC = %.3f\n', nanmean(icc_vals));
fprintf('Median ICC = %.3f\n', nanmedian(icc_vals));
fprintf('Range = [%.3f, %.3f]\n', nanmin(icc_vals), nanmax(icc_vals));

% 可视化 (直方图)
figure;
histogram(icc_vals, 20);
xlabel('ICC(2,1)');
ylabel('Number of ROIs');
title('Distribution of ROI-wise ICC values');
grid on;

% 如果你有脑图 parcellation，可以把 icc_vals 映射回脑表面显示
% （这里就留个占位，需用你已有的绘图函数，比如 plot_cortical 或 FreeSurfer fsaverage）






%% Bland–Altman cbf
% 假设你的数据是两个向量:
test = CBF1_ValPC1_bil(CBFsnr_mean_bil>0);     %[nROI x 1]  (test session 的 CCSI)
retest = CBF2_ValPC1_bil(CBFsnr_mean_bil>0);   %[nROI x 1] (retest session 的 CCSI)
% test = mean(Data1_Val_bil(CBFsnr_mean_bil>0),2);     %[nROI x 1]  (test session 的 CCSI)
% retest = mean(Data2_Val_bil(CBFsnr_mean_bil>0),2);   %[nROI x 1] (retest session 的 CCSI)

% 1. 计算均值和差值
meanVals = (test + retest) / 2;   % 横轴
diffVals = test - retest;         % 纵轴

% 2. 计算 bias 和一致性限
bias = mean(diffVals);            
sd_diff = std(diffVals);          
loa_upper = bias + 1.96 * sd_diff; % 上一致性限
loa_lower = bias - 1.96 * sd_diff; % 下一致性限

% 3. 作图
figure; hold on;
scatter(meanVals, diffVals, 40, 'k', 'filled'); % 散点图
yline(bias, 'r-', 'LineWidth', 2, 'Label', 'Bias'); 
yline(loa_upper, 'b--', 'LineWidth', 2, 'Label', '+1.96 SD');
yline(loa_lower, 'b--', 'LineWidth', 2, 'Label', '-1.96 SD');
xlabel('Mean of Test and Retest CBF');
ylabel('Difference (Test - Retest)');
title('Bland-Altman Plot for Test–Retest CBF');
grid on;
set(gca, 'FontSize', 12);

% 4. 输出结果
fprintf('Bias = %.4f\n', bias);
fprintf('Upper LoA = %.4f\n', loa_upper);
fprintf('Lower LoA = %.4f\n', loa_lower);



%% roi-wise icc cbf
% 假设你的数据矩阵如下：
% test_data  = [nSubjects x nROIs]   % test session 的 CCSI
% retest_data = [nSubjects x nROIs]  % retest session 的 CCSI
% nSubjects = 14, nROIs = 360 (HCP-MMP1 parcellation)
test_data = Data1_Val_bil(CBFsnr_mean_bil>0,:)';     %[nROI x 1]  (test session 的 CCSI)
retest_data = Data2_Val_bil(CBFsnr_mean_bil>0,:)';   %[nROI x 1] (retest session 的 CCSI)

nSubjects = size(test_data,1);
nROIs = size(test_data,2);

icc_vals = nan(nROIs,1);

% 逐 ROI 计算 ICC(2,1) - 双向随机效应，单次测量，绝对一致性
for roi = 1:nROIs
    Y = [test_data(:,roi), retest_data(:,roi)];  % [nSubjects x 2]
    
    % 去掉 NaN
    Y(any(isnan(Y),2),:) = [];
    
    % 计算方差分解
    n = size(Y,1); % 有效样本数
    k = size(Y,2); % 重复次数 (这里是2: test & retest)
    mean_per_subject = mean(Y,2);
    mean_per_rater   = mean(Y,1);
    grand_mean = mean(Y(:));
    
    % 平方和
    SS_between = k * sum((mean_per_subject - grand_mean).^2);
    SS_within  = sum(sum((Y - mean_per_subject).^2));
    SS_rater   = n * sum((mean_per_rater - grand_mean).^2);
    
    % 均方
    MS_between = SS_between / (n-1);
    MS_within  = SS_within / ((n-1)*(k-1));
    MS_rater   = SS_rater / (k-1);
    
    % ICC(2,1) 公式
    icc = (MS_between - MS_within) / ...
          (MS_between + (k-1)*MS_within + (k/n)*(MS_rater - MS_within));
    
    icc_vals(roi) = icc;
end

% 输出结果
fprintf('ROI-wise ICC summary:\n');
fprintf('Mean ICC = %.3f\n', nanmean(icc_vals));
fprintf('Median ICC = %.3f\n', nanmedian(icc_vals));
fprintf('Range = [%.3f, %.3f]\n', nanmin(icc_vals), nanmax(icc_vals));

% 可视化 (直方图)
figure;
histogram(icc_vals, 20);
xlabel('ICC(2,1)');
ylabel('Number of ROIs');
title('Distribution of ROI-wise ICC values');
grid on;

% 如果你有脑图 parcellation，可以把 icc_vals 映射回脑表面显示
% （这里就留个占位，需用你已有的绘图函数，比如 plot_cortical 或 FreeSurfer fsaverage）

%% batch write to 1D for compute spin-null by neuromaps
% bil
OutNames = {'CBF_Mean_PScs', 'CBF_PC1_Val',...
            'CBF_PC1_Val_Layer1',...
            'CBF_PC1_Val_Layer2',...
            'CBF_PC1_Val_Layer3',...
            'CBF_PC1_Val_Layer4',...
            'CBF_PC1_Val_Layer5',...
            'CBF_PC1_Val_Layer6'};  % there layers 1-6 actually are L6-L1
OutDatas = {CBF_meanPScs_bil, CBF_ValPC1_bil,...
            CBF_ValPC1_Layer_bil{1},...
            CBF_ValPC1_Layer_bil{2},...
            CBF_ValPC1_Layer_bil{3},...
            CBF_ValPC1_Layer_bil{4},...
            CBF_ValPC1_Layer_bil{5},...
            CBF_ValPC1_Layer_bil{6},...
            CBFsnr_mean_bil, CNRcm_g};
UniOrBil = 2;  % 1:uni  2:bil
for i = 1:numel(OutNames)
    if UniOrBil == 1
        OutData = [OutDatas{i}; OutDatas{i}];
        suffix = '_uni.1D';
    else
        OutData = OutDatas{i};
        suffix = '_bil.1D'
    end
    % do it
    OutDir = '../Data/group';
    save([OutDir '/' OutNames{i} suffix],'OutData','-ascii');
end


% uni
OutNames = {'CBF_Mean_PScs', 'CBF_PC1_Val',...
            'CBF_PC1_Val_Layer1',...
            'CBF_PC1_Val_Layer2',...
            'CBF_PC1_Val_Layer3',...
            'CBF_PC1_Val_Layer4',...
            'CBF_PC1_Val_Layer5',...
            'CBF_PC1_Val_Layer6'};  % there layers 1-6 actually are L6-L1
OutDatas = {CBF_meanPScs_uni, CBF_ValPC1_uni,...
            CBF_ValPC1_Layer_uni{1},...
            CBF_ValPC1_Layer_uni{2},...
            CBF_ValPC1_Layer_uni{3},...
            CBF_ValPC1_Layer_uni{4},...
            CBF_ValPC1_Layer_uni{5},...
            CBF_ValPC1_Layer_uni{6}};
UniOrBil = 1;  % 1:uni  2:bil
for i = 1:numel(OutNames)
    if UniOrBil == 1
        OutData = [OutDatas{i}; OutDatas{i}];
        suffix = '_uni.1D';
    else
        OutData = OutDatas{i};
        suffix = '_bil.1D';
    end
    % do it
    OutDir = '../Data/group';
    save([OutDir '/' OutNames{i} suffix],'OutData','-ascii');
end


%% write to map plot by R
OutNames = {'CBF_Mean_PScs_bil.csv', 'CBF_PC1_Val_bil.csv',...
            'Mito_CI_bil.csv', 'Mito_CII_bil.csv',...
            'Mito_CIV_bil.csv', 'Mito_MitoD_bil.csv',...
            'BB_G1_bil.csv', 'CBF_PC1_Val_L6_bil.csv',...
            'CBF_PC1_Val_L5_bil.csv', 'CBF_PC1_Val_L4_bil.csv',...
            'CBF_PC1_Val_L3_bil.csv', 'CBF_PC1_Val_L2_bil.csv',...
            'CBF_PC1_Val_L1_bil.csv', 'intersubjvar.csv',...
            'myelinmap.csv', 'thickness.csv',...
            'Economo-Koskinas-cytoarchitectonics.csv', 'StructualConnection_G1.csv',...
            'Yeo2011-FunctionNetwork.csv', 'CNR_G1.csv',...
            'SNRmean.csv'};
OutDatas = {CBF_meanPScs_bil, CBF_ValPC1_bil,...
            zMitoData(:,1), zMitoData(:,2),...
            zMitoData(:,3), zMitoData(:,4),...
            BB_gradient, CBF_ValPC1_Layer_bil{1},...
            CBF_ValPC1_Layer_bil{2}, CBF_ValPC1_Layer_bil{3},...
            CBF_ValPC1_Layer_bil{4}, CBF_ValPC1_Layer_bil{5},...
            CBF_ValPC1_Layer_bil{6}, all_other_maps.maps{48},...
            all_other_maps.maps{26}, all_other_maps.maps{27},...
            EKc, scm_g(:,1),...
            Yeo_FN, CNRcm_g(:,1),...
            CBFsnr_mean_bil};
Datalabel = [1:360]';

for i = 1:numel(OutNames)
    tmp = subfun_write_to_MapPlot_byR(OutDatas{i},Datalabel,atlas_info,[Group_dir OutNames{i}]);
end




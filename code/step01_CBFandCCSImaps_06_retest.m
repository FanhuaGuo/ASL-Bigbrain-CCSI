clc
clear all
close all
set(0,'defaultfigurecolor',[1 1 1]);% set figure background color
addpath(genpath('/Users/guofanhua/Desktop/gfh/tools/MatlabTools/BrainSpace-latest/matlab'))

%% ======================= step 1 for disposal all data

%% input data
step01_CBFandCCSImaps_00_retest_init;
step01_CBFandCCSImaps_00_init;
step01_CBFandCCSImaps_00_ComputeBasicResults;

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



%% regional 14-retest CCSI with SNR>4 threshold (Fig S5)
[~,p] = ttest(Data2_PScs_uni');
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
data = Data2_PScs_uni(CBFsnr_mean_uni>4,:)';
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
title('ROI-based Bland-Altman Plot for 14 Test–Retest CCSI');
grid on;
set(gca, 'FontSize', 12);

% 4. 输出结果
fprintf('Bias = %.4f\n', bias);
fprintf('Upper LoA = %.4f\n', loa_upper);
fprintf('Lower LoA = %.4f\n', loa_lower);
fprintf('(2 * 1.96 * sd_diff)/(Meanmax - Meanmin) = %.4f\n', (loa_upper-loa_lower)/(max(meanVals)-min(meanVals)));




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
xlabel('Mean of Test and Retest CBF-score');
ylabel('Difference (Test - Retest)');
title('ROI-based Bland-Altman Plot for 14 Test–Retest CBF-score');
grid on;
set(gca, 'FontSize', 12);

% 4. 输出结果
fprintf('Bias = %.4f\n', bias);
fprintf('Upper LoA = %.4f\n', loa_upper);
fprintf('Lower LoA = %.4f\n', loa_lower);
fprintf('(2 * 1.96 * sd_diff)/(Meanmax - Meanmin) = %.4f\n', (loa_upper-loa_lower)/(max(meanVals)-min(meanVals)));



%% batch write to 1D for compute spin-null by neuromaps
% bil
OutNames = {'retest_Mean_PScs'};  % there layers 1-6 actually are L6-L1
OutDatas = {CBF2_meanPScs_bil};
UniOrBil = 2;  % 1:uni  2:bil
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


% uni
OutNames = {'retest_Mean_PScs'};  % there layers 1-6 actually are L6-L1
OutDatas = {CBF2_meanPScs_uni};
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
OutNames = {'CBF_Mean_PScs2_qc_bil.csv'};
SNRmean = CBFsnr_mean_bil;
OutDatas = {CBF2_meanPScs_bil(SNRmean>4)};
Datalabel = [1:360]';
Datalabel(SNRmean<=4) = [];

for i = 1:numel(OutNames)
    tmp = subfun_write_to_MapPlot_byR(OutDatas{i},Datalabel,atlas_info,[Group_dir OutNames{i}]);
end





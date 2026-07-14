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


%% CBF PC2/3
[tmp, explained] = subfun_PCAtoPC(zscore(Data_Val_bil));
CBF_ValPC1_bil = tmp(:,1);
CBF_ValPC2_bil = tmp(:,2);
CBF_ValPC3_bil = tmp(:,3);


[tmp, explained_uni] = subfun_PCAtoPC(zscore(Data_Val_uni));
CBF_ValPC1_uni = tmp(:,1);
CBF_ValPC2_uni = tmp(:,2);
CBF_ValPC3_uni = tmp(:,3);

%% lobe map
tmp = atlas_info.Lobe;
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
atlas_lobeindex(181:end) = atlas_lobeindex(1:180);

%% check atlas label
check = [atlas_info.regionLongName(1:180) atlas_info.regionLongName(181:360) atlas_info.Lobe(1:180) atlas_info.Lobe(181:360)];


%% plot BB/CBF LP whole-brain
figure('Color',[1 1 1],'Position',[0 0 500 500],'Units','pixels');
hold on
tmpPlot = zscore(BB_LaminarProfile(:,1:end-1)')';
for i = 1:size(tmpPlot,1)
    plot(tmpPlot(i,:),'-','linewidth',0.5);
end

figure('Color',[1 1 1],'Position',[0 0 400 300],'Units','pixels');
hold on
tmpPlot = zeros(360,12);
for i = 1:numel(Data_scan)
    tmpPlot = tmpPlot + zscore(Data_scan{i}(:,1:end)')';
end
tmpPlot = tmpPlot/numel(Data_scan);
for i = 1:size(tmpPlot,1)
    plot(tmpPlot(i,:),'-','linewidth',0.5);
end
set(gca,'ylim',[-2.5 2], 'yTick',[])
set(gca,'xlim',[0 13], 'xTick',[])
xlabel('cortical depth')
ylabel('CBF')
box on
% ylim([-2.5 2]); yTick([]);


%% write to SourceData
% Fig 2ade
outSD = [atlas_info.regionName atlas_info.regionLongName];
tmp = table(CBF_Valmean_bil, CBF_ValPC1_bil, CBF_ValPC1_Layer_bil{1}, CBF_ValPC1_Layer_bil{2}, CBF_ValPC1_Layer_bil{3},...
            'VariableNames',{'CBF mean', 'CBF PC1', 'CBF PC1 L1', 'CBF PC1 L2', 'CBF PC1 L3'});
outSD = [outSD tmp];


% Fig 3ce
outSD = [atlas_info.regionName atlas_info.regionLongName];
tmp = table(CBF_meanPScs_bil_correct, T1vBB_meanPScs_bil_correct,...
            'VariableNames',{'CCSI', 'CCSIcontrol'});
outSD = [outSD tmp];
outSD.CCSI(CBFsnr_mean_bil<=4) = nan;
outSD.CCSIcontrol(CBFsnr_mean_bil<=4) = nan;


% Extended Data Fig 3adf
outSD = [atlas_info.regionName atlas_info.regionLongName];
tmp = table(CBF_meanPScs_bil, CBF_meanPScs_bil_correct, CBFsnr_mean_bil,...
            'VariableNames',{'origCCSI', 'correctedCCSI', 'SNR'});
outSD = [outSD tmp];




%% batch write to 1D for compute spin-null by neuromaps
suffix0 = ['_L' num2str(num_layer)];

% bil
OutNames = {'CBF_Mean_PScs', 'CBF_PC1_Val', 'CBF_mean_Val'}; 
OutDatas = {CBF_meanPScs_bil_correct, CBF_ValPC1_bil, CBF_Valmean_bil};
% OutDatas = {CBF_meanPScs_bil, CBF_ValPC1_bil, CBF_Valmean_bil};
if num_layer == 12
    OutNames = {OutNames{:}, ...
                'T1vCBF_Mean_PScs', 'T1vBB_Mean_PScs',...
                'T1wCBF_Mean_PScs', 'T1wBB_Mean_PScs',...
                'CBF_Mean_PScs_uncorr'};
    OutDatas = {OutDatas{:}, ...
                T1vCBF_meanPScs_bil_correct, T1vBB_meanPScs_bil_correct,... 
                T1wCBF_meanPScs_bil_correct, T1wBB_meanPScs_bil_correct,...
                CBF_meanPScs_bil};
%     OutDatas = {OutDatas{:}, ...
%                 T1vCBF_meanPScs_bil, T1vBB_meanPScs_bil,... 
%                 T1wCBF_meanPScs_bil, T1wBB_meanPScs_bil};
end
UniOrBil = 2;  % 1:uni  2:bil
for i = 1%:numel(OutNames)
    if UniOrBil == 1
        OutData = [OutDatas{i}; OutDatas{i}];
        suffix = '_uni.1D';
    else
        OutData = OutDatas{i};
        suffix = '_bil.1D';
    end
    % do it
    OutDir = '../Data/group';
    save([OutDir '/' OutNames{i} suffix0 suffix],'OutData','-ascii');
end


% uni
OutNames = {'CBF_Mean_PScs', 'CBF_PC1_Val', 'CBF_mean_Val'};  
OutDatas = {CBF_meanPScs_uni, CBF_ValPC1_uni, CBF_Valmean_uni};
if num_layer == 12
    OutNames = {OutNames{:}, ...
                'T1vCBF_Mean_PScs', 'T1vBB_Mean_PScs',...
                'T1wCBF_Mean_PScs', 'T1wBB_Mean_PScs'};
    OutDatas = {OutDatas{:}, ...
                T1vCBF_meanPScs_uni, T1vBB_meanPScs_uni,...
                T1wCBF_meanPScs_uni, T1wBB_meanPScs_uni};
end
UniOrBil = 1;  % 1:uni  2:bil
for i = 1%:numel(OutNames)
    if UniOrBil == 1
        OutData = [OutDatas{i}; OutDatas{i}];
        suffix = '_uni.1D';
    else
        OutData = OutDatas{i};
        suffix = '_bil.1D';
    end
    % do it
    OutDir = '../Data/group';
    save([OutDir '/' OutNames{i} suffix0 suffix],'OutData','-ascii');
end


%% write to map plot by R MRI maps
suffix0 = ['_L' num2str(num_layer)];
OutNames = {'CBF_Mean_PScs_bil', 'CBF_PC1_Val_bil',...
            'CBF_mean_Val_bil', 'SNRmean',...
            'CBF_PC1_Val_deep_bil', 'CBF_PC1_Val_middle_bil', ...
            'CBF_PC1_Val_superficial_bil', 'CBF_Mean_PScs_corr_bil'};
OutDatas = {CBF_meanPScs_bil, CBF_ValPC1_bil,...
            CBF_Valmean_bil, CBFsnr_mean_bil,...
            CBF_ValPC1_Layer_bil{1}, CBF_ValPC1_Layer_bil{2},...
            CBF_ValPC1_Layer_bil{3}, CBF_meanPScs_bil_correct};
if num_layer == 12
%     OutNames = {OutNames{:}, ...
%                 'T1vCBF_Mean_PScs_bil', 'T1vBB_Mean_PScs_bil',...
%                 'T1wCBF_Mean_PScs_bil', 'T1wBB_Mean_PScs_bil'};
%     OutDatas = {OutDatas{:}, ...
%                 T1vCBF_meanPScs_bil, T1vBB_meanPScs_bil,...
%                 T1wCBF_meanPScs_bil, T1wBB_meanPScs_bil};
    OutNames = {OutNames{:}, ...
                'T1vCBF_Mean_PScs_corr_bil', 'T1vBB_Mean_PScs_corr_bil',...
                'T1wCBF_Mean_PScs_corr_bil', 'T1wBB_Mean_PScs_corr_bil'};
    OutDatas = {OutDatas{:}, ...
                T1vCBF_meanPScs_bil_correct, T1vBB_meanPScs_bil_correct,... 
                T1wCBF_meanPScs_bil_correct, T1wBB_meanPScs_bil_correct};
end

Datalabel = [1:360]';

for i = 1:numel(OutNames)
    tmp = subfun_write_to_MapPlot_byR(OutDatas{i},Datalabel,atlas_info,[Group_dir OutNames{i} suffix0 '.csv']);
end




%% write to map plot by R other maps
OutNames = {'Mito_CI_bil.csv', 'Mito_CII_bil.csv',...
            'Mito_CIV_bil.csv', 'Mito_MitoD_bil.csv',...
            'Mito_MRC_bil.csv', 'Mito_TRC_bil.csv',...
            'BB_G1_bil.csv', 'intersubjvar.csv',...
            'myelinmap.csv', 'thickness.csv',...
            'Economo-Koskinas-cytoarchitectonics.csv', 'StructualConnection_G1.csv',...
            'Yeo2011-FunctionNetwork.csv', 'FunctionGradient1.csv',...
            'lobe.csv', 'FunctionGradient1_enigma.csv',...
            'qT1Gradient1.csv', 'T1wGradient1.csv',...
            'PET_CBF.csv'};
OutDatas = {zMitoData(:,1), zMitoData(:,2),...
            zMitoData(:,3), zMitoData(:,4),...
            zMitoData(:,5), zMitoData(:,6),...
            BB_gradient, all_other_maps.maps{48},...
            all_other_maps.maps{26}, all_other_maps.maps{27},...
            EKc, scm_g(:,1),...
            Yeo_FN, all_other_maps.maps{38},...
            atlas_lobeindex, -fcm_g(:,1),...
            -qT1_g(:,1), T1w_g(:,1),...
            all_other_maps.maps{54}};
Datalabel = [1:360]';

for i = 1:numel(OutNames)
    tmp = subfun_write_to_MapPlot_byR(OutDatas{i},Datalabel,atlas_info,[Group_dir OutNames{i}]);
end




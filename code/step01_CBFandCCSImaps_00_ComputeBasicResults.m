
%% ======================= step 01-00 for disposal all data
% compute some basic maps

%% compute cosine similarity between CBF and Cell-body staining intensity (CCSI)
Data_PScs_bil = subfun_ProfileSimilarity(Data_scan_bil,Data_BBLP_bil);
Data_PScs_uni = subfun_ProfileSimilarity(Data_scan_uni,Data_BBLP_uni);

Data_PScs_bil = pi-acos(Data_PScs_bil) -pi/2;
Data_PScs_uni = pi-acos(Data_PScs_uni) -pi/2;


if num_layer == 12
    Data_PScs_T1wCBF_bil = subfun_ProfileSimilarity_mri2mri(Data_scan_bil,Data_T1w_bil);
    Data_PScs_T1wCBF_bil = pi-acos(Data_PScs_T1wCBF_bil) -pi/2;
    Data_PScs_T1vCBF_bil = subfun_ProfileSimilarity_mri2mri(Data_scan_bil,Data_T1v_bil);
    Data_PScs_T1vCBF_bil = pi-acos(Data_PScs_T1vCBF_bil) -pi/2;
    Data_PScs_T1wCBF_uni = subfun_ProfileSimilarity_mri2mri(Data_scan_uni,Data_T1w_uni);
    Data_PScs_T1wCBF_uni = pi-acos(Data_PScs_T1wCBF_uni) -pi/2;
    Data_PScs_T1vCBF_uni = subfun_ProfileSimilarity_mri2mri(Data_scan_uni,Data_T1v_uni);
    Data_PScs_T1vCBF_uni = pi-acos(Data_PScs_T1vCBF_uni) -pi/2;
    
    Data_PScs_T1wBB_bil = subfun_ProfileSimilarity(Data_T1w_bil,Data_BBLP_bil);
    Data_PScs_T1wBB_bil = pi-acos(Data_PScs_T1wBB_bil) -pi/2;
    Data_PScs_T1vBB_bil = subfun_ProfileSimilarity(Data_T1v_bil,Data_BBLP_bil);
    Data_PScs_T1vBB_bil = pi-acos(Data_PScs_T1vBB_bil) -pi/2;
    Data_PScs_T1wBB_uni = subfun_ProfileSimilarity(Data_T1w_uni,Data_BBLP_uni);
    Data_PScs_T1wBB_uni = pi-acos(Data_PScs_T1wBB_uni) -pi/2;
    Data_PScs_T1vBB_uni = subfun_ProfileSimilarity(Data_T1v_uni,Data_BBLP_uni);
    Data_PScs_T1vBB_uni = pi-acos(Data_PScs_T1vBB_uni) -pi/2;
end


% correct left-right bias (Global Mean Offset Correction)
Data_PScs_bil_correct = subfun_GlobalMeanOffsetCorrection(Data_PScs_bil);
if num_layer == 12
    Data_PScs_T1wCBF_bil_correct = subfun_GlobalMeanOffsetCorrection(Data_PScs_T1wCBF_bil);
    Data_PScs_T1vCBF_bil_correct = subfun_GlobalMeanOffsetCorrection(Data_PScs_T1vCBF_bil);
    Data_PScs_T1wBB_bil_correct = subfun_GlobalMeanOffsetCorrection(Data_PScs_T1wBB_bil);
    Data_PScs_T1vBB_bil_correct = subfun_GlobalMeanOffsetCorrection(Data_PScs_T1vBB_bil);
end



%% mean or PCA combine all subjects into 1 map
% cbf
CBF_Valmean_bil = mean(Data_Val_bil,2);
CBF_Valmean_uni = mean(Data_Val_uni,2);

tmp = subfun_PCAtoPC(zscore(Data_Val_uni));     % whole-brain PCA
CBF_ValPC1_uni = tmp(:,1);

[tmp, explained] = subfun_PCAtoPC(zscore(Data_Val_bil));  % whole-brain PCA
CBF_ValPC1_bil = tmp(:,1);

% [tmp, explained] = subfun_PCAtoPC(zscore(Data_Val_bil(1:end/2,:)));  % half-brain PCA
% CBF_ValPC1_bil = tmp(:,1);
% [tmp, explained] = subfun_PCAtoPC(zscore(Data_Val_bil(end/2+1:end,:))); 
% CBF_ValPC1_bil = [CBF_ValPC1_bil; tmp(:,1)];


% laminar cbf -- 3 layers
CBF_ValPC1_Layer_uni = {};
CBF_ValPC1_Layer_bil = {};

for i = 1:3
    tmp = subfun_PCAtoPC(zscore(Data_ValLayer_uni{i}));
    CBF_ValPC1_Layer_uni{i} = tmp(:,1);
    tmp = subfun_PCAtoPC(zscore(Data_ValLayer_bil{i}));
    CBF_ValPC1_Layer_bil{i} = tmp(:,1);
end



% ccsi
CBF_meanPScs_bil = mean(Data_PScs_bil,2);
CBF_meanPScs_uni = mean(Data_PScs_uni,2);

CBF_meanPScs_bil_correct = mean(Data_PScs_bil_correct,2);

CBFsnr_mean_bil = mean(Data_SNR,2);
CBFsnr_mean_uni = subfun_merge_bil2uni(CBFsnr_mean_bil);


% T1
if num_layer == 12
    T1wCBF_meanPScs_bil = mean(Data_PScs_T1wCBF_bil,2);
    T1vCBF_meanPScs_bil = mean(Data_PScs_T1vCBF_bil,2);
    T1wCBF_meanPScs_uni = mean(Data_PScs_T1wCBF_uni,2);
    T1vCBF_meanPScs_uni = mean(Data_PScs_T1vCBF_uni,2);
    T1wCBF_meanPScs_bil_correct = mean(Data_PScs_T1wCBF_bil_correct,2);
    T1vCBF_meanPScs_bil_correct = mean(Data_PScs_T1vCBF_bil_correct,2);
    
    T1wBB_meanPScs_bil = mean(Data_PScs_T1wBB_bil,2);
    T1vBB_meanPScs_bil = mean(Data_PScs_T1vBB_bil,2);
    T1wBB_meanPScs_uni = mean(Data_PScs_T1wBB_uni,2);
    T1vBB_meanPScs_uni = mean(Data_PScs_T1vBB_uni,2);
    T1wBB_meanPScs_bil_correct = mean(Data_PScs_T1wBB_bil_correct,2);
    T1vBB_meanPScs_bil_correct = mean(Data_PScs_T1vBB_bil_correct,2);
end



%% gradient compute--structual, functional gradient and CNR matrix gradient
fcm = load([Refer_dir 'enigma/funcMatrix_ctx_glasser_360.csv']);  % data from enigmatools
gm = GradientMaps('kernel','cs','approach','dm','n_components',10,'alignment','none'); % none procrustes
gm = gm.fit(fcm);
fcm_g = gm.gradients{1};  % did not use, I directly used the fcgradient1-10 in neuromaps 

scm = load([Refer_dir 'enigma/strucMatrix_ctx_glasser_360.csv']);  % data from enigmatools
gm = GradientMaps('kernel','cs','approach','dm','n_components',10,'alignment','none');
gm = gm.fit(scm);
scm_g = gm.gradients{1};


% CNR gradient
CNRcm = subfun_merge_CNRm2CNRcm(Data_CNRm);
CNRcm(1:size(CNRcm,1)+1:end) = 0; 
CNRcm(isnan(CNRcm)) = 0;
gm = GradientMaps('kernel','cs','approach','dm','n_components',10,'alignment','none');
gm = gm.fit(CNRcm);
CNRcm_g = gm.gradients{1};


%% gradient compute--cbf 
tmpMeanm = zeros(360,360);
for i = 1:numel(Data_scan_bil)
%     tmp = corr(Data_scan_bil{i}', 'Type', 'Pearson', 'Rows', 'pairwise');
%     tmpMeanm = tmpMeanm + subfun_fisher_z(tmp);
    
    tmpMeanm = tmpMeanm + subfun_computeMPC(Data_scan_bil{i});
end
cbfm = tmpMeanm/numel(Data_scan_bil);
% cbfm = subfun_fisher_z_inverse(cbfm);
% cbfm(1:361:end) = 1;


gm = GradientMaps('kernel','cs','approach','dm','n_components',10,'alignment','none');
gm = gm.fit(cbfm);
cbfm_g = gm.gradients{1};  % did not use, I directly used the fcgradient1-10 in neuromaps 



% corr(fcm_g,scm_g(:,1:end-1))



%% gradient compute--T1w
tmpMeanm = zeros(360,360);
for i = 1:numel(Data_T1w_bil)
%     tmp = corr(Data_T1w_bil{i}', 'Type', 'Pearson', 'Rows', 'pairwise');
%     tmpMeanm = tmpMeanm + subfun_fisher_z(tmp);
    
    tmpMeanm = tmpMeanm + subfun_computeMPC(Data_T1w_bil{i});
end
T1wm = tmpMeanm/numel(Data_T1w_bil);
% T1wm = subfun_fisher_z_inverse(T1wm);
% T1wm(1:361:end) = 1;


gm = GradientMaps('kernel','cs','approach','dm','n_components',10,'alignment','none');
gm = gm.fit(T1wm);
T1w_g = gm.gradients{1};  % did not use, I directly used the fcgradient1-10 in neuromaps 


% corr(-fcm_g(:,1),qT1_g(:,1:end-1),'Type','Spearman')
% corr(-fcm_g(:,1),scm_g(:,1:end-1))
% corr(BB_gradient,qT1_g(:,1:end-1))
% corr(BB_gradient,T1w_g(:,1:end-1))
% corr(BB_gradient,-fcm_g(:,1),'Type','Spearman') % Pearson Spearman
% corr(-fcm_g(:,1),all_other_maps.maps{38})


%% MPC gradient compute--T1v
tmpMeanm = zeros(360,360);
for i = 1:numel(Data_T1v_bil)
%     tmp = corr(Data_T1v_bil{i}', 'Type', 'Pearson', 'Rows', 'pairwise');
%     tmpMeanm = tmpMeanm + subfun_fisher_z(tmp);
    
    tmpMeanm = tmpMeanm + subfun_computeMPC(Data_T1v_bil{i});
end
qT1m = tmpMeanm/numel(Data_T1v_bil);
% qT1m = subfun_fisher_z_inverse(tmpMeanm);
% qT1m(1:361:end) = 1;


gm = GradientMaps('kernel','cs','approach','dm','n_components',10,'alignment','none');
gm = gm.fit(qT1m);
qT1_g = gm.gradients{1};  % did not use, I directly used the fcgradient1-10 in neuromaps 

% %%
% % 1. 获取排序索引 (按 G1 的值升序排列)
% % 注意：如果你使用的是 neuromaps 的 G1，确保其长度与矩阵维度(360)一致
% [~, sort_idx] = sort(qT1_g(:,1), 'ascend');
% 
% % 2. 准备要可视化的矩阵
% % 论文中展示的 Figure B 通常是 Affinity Matrix (关联矩阵)
% % 如果你使用了 BrainSpace 的 GradientMaps，关联矩阵存在 gm.af{1} 中
% matrix_to_show = qT1m; 
% 
% 
% % 3. 对矩阵的行和列同时进行重排
% sorted_matrix = matrix_to_show(sort_idx, sort_idx);
% 
% % 4. 绘图可视化
% % figure('Color', 'w');
% % imagesc(sorted_matrix);
% % colormap(jet); % 或者使用类似论文的 'magma' 颜色图
% % colorbar;
% % axis square;
% % title('G1 Ordered Matrix');
% % xlabel('Nodes (Sorted by G1)');
% % ylabel('Nodes (Sorted by G1)');
% 
% % 建议：为了让对比更明显，可以调整显示的数值范围
% % caxis([0, 1]); % 如果是归一化角度矩阵
% 
% 




%% ======================= step 01-00 for disposal all data
% compute some basic maps

%% compute cosine similarity between CBF and Cell-body staining intensity (CCSI)
Data_PScs_bil = subfun_ProfileSimilarity(Data_scan_bil,Data_BBLP_bil);
Data_PScs_uni = subfun_ProfileSimilarity(Data_scan_uni,Data_BBLP_uni);

Data_PScs_bil = pi-acos(Data_PScs_bil) -pi/2;
Data_PScs_uni = pi-acos(Data_PScs_uni) -pi/2;


%% mean or PCA combine all subjects into 1 map
tmp = subfun_PCAtoPC(zscore(Data_Val_uni));
CBF_ValPC1_uni = tmp(:,1);
[tmp, explained] = subfun_PCAtoPC(zscore(Data_Val_bil));
CBF_ValPC1_bil = tmp(:,1);

CBF_ValPC1_Layer_uni = {};
CBF_ValPC1_Layer_bil = {};
for i = 1:num_layer
    tmp = subfun_PCAtoPC(zscore(Data_ValLayer_uni{i}));
    CBF_ValPC1_Layer_uni{i} = tmp(:,1);
    tmp = subfun_PCAtoPC(zscore(Data_ValLayer_bil{i}));
    CBF_ValPC1_Layer_bil{i} = tmp(:,1);
end

CBF_meanPScs_bil = mean(Data_PScs_bil,2);
CBF_meanPScs_uni = mean(Data_PScs_uni,2);

CBFsnr_mean_bil = mean(Data_SNR,2);
CBFsnr_mean_uni = subfun_merge_bil2uni(CBFsnr_mean_bil);




%% gradient compute--structual, functional gradient and CNR matrix gradient
fcm = load([Refer_dir 'enigma/funcMatrix_ctx_glasser_360.csv']);  % data from enigmatools
gm = GradientMaps('kernel','cs','approach','dm','n_components',10,'alignment','none');
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




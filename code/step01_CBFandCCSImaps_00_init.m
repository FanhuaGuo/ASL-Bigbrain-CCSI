
%% ======================= step 01-00 for disposal all data
% input all the basic data

%% set parameters
% input parameters
subjs = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10',...
         'S11','S12','S13','S14','S15','S16','S17','S18','S19','S20',...
         'S21','S22','S23','S24','S25','S26','S27','S28','S29','S30',...
         'S31','S32','S33','S34','S35','S36','S37','S38','S39','S40',...
         'S41','S42','S43','S44','S45','S46','S47','S48','S49','S50'};
% which_subjs = [1:9 13 15:22];
which_subjs = [1:9 13 15:22 32 36:46];
% which_subjs = [1:6 8 9 13 15:22 32 36:46];
num_layer = 12;     % 6 10 12 15 20
atlas = 'glasser360std141';     % glasser360  glasser360std141
ifloadGO = 1;
signal = 'CBF';
FileName_Mito = ['Mito_glasser360_bil.1D'];
scan_suffix = '_s1';

% basic parameters
Refer_dir = '../../reference/';
Bigbrain_dir = [Refer_dir 'BigBrain/'];
Data_dir = '../Data/';
Group_dir = [Data_dir 'group/'];
CellType_dir = [Refer_dir 'CellType/'];
FileName_BigBrain = ['wholeBB_glasser360_' num2str(num_layer) 'layers_LaminarProfile.mat'];



%% read scan and BigBrain data
% read scan data
Data_T1w = {};
Data_T1v = {};
tmpT1v = zeros(360,num_layer);
Data_scan = {};
Data_Val = [];
Data_ValLayer = {}; Data_ValLayer{num_layer} = [];
Data_SNR = [];
Data_CNRm = {};
Data_ROIvol = [];
for i = which_subjs
    subj = subjs{i};
    tmp = load([Data_dir subj '/LaminarProfile_' num2str(num_layer) 'points_' atlas '_' signal scan_suffix '.mat']);
    Data_scan = {Data_scan{:}, tmp.LaminarProfile};
    Data_Val = [Data_Val tmp.MRIsignal];
    for k = 1:3
        tmplayers = round((k-1)/3*num_layer)+1 : round(k/3*num_layer);
        Data_ValLayer{k} = [Data_ValLayer{k} mean(tmp.LaminarProfile(:,tmplayers),2)];
    end
    Data_SNR = [Data_SNR tmp.SNR];
    Data_CNRm{size(Data_SNR,2)} = tmp.CNRm;
    Data_ROIvol = [Data_ROIvol sum(tmp.LaminarProfile_voxelsNum,2)*0.5^3];
    
    if num_layer == 12
        tmp = load([Data_dir subj '/LaminarProfile_' num2str(num_layer) 'points_' atlas '_T1wv.mat']);
        Data_T1w = {Data_T1w{:}, tmp.LaminarProfileT1w};
        Data_T1v = {Data_T1v{:}, tmp.LaminarProfileT1v};
        
        if i~=2
            tmpT1v = tmpT1v+tmp.LaminarProfileT1v;
        end
    end
end

if num_layer == 12
    Data_T1v{2} = tmpT1v/(numel(which_subjs)-1);
    Data_T1w_bil = Data_T1w';
    Data_T1v_bil = Data_T1v';
    Data_T1w_uni = {};
    Data_T1v_uni = {};
    for i = 1:numel(Data_T1w_bil)
        Data_T1w_uni{i} = subfun_merge_bil2uni(Data_T1w_bil{i}')';
        Data_T1v_uni{i} = subfun_merge_bil2uni(Data_T1v_bil{i}')';
    end
end
Data_Val_bil = Data_Val;
Data_Val_uni = subfun_merge_bil2uni(Data_Val_bil')';
Data_ValLayer_bil = Data_ValLayer;
Data_ValLayer_uni = {};
for i = 1:3 %num_layer
    Data_ValLayer_uni{i} = subfun_merge_bil2uni(Data_ValLayer_bil{i}')';
end

Data_scan_bil = Data_scan';
Data_scan_uni = {};
for i = 1:numel(Data_scan_bil)
    Data_scan_uni{i} = subfun_merge_bil2uni(Data_scan_bil{i}')';
end



%% read bigbrain data
load([Bigbrain_dir FileName_BigBrain]);
Data_BBLP_bil = BB_LaminarProfile;
Data_BBLP_uni = subfun_merge_bil2uni(Data_BBLP_bil')';



%% read atlas info
atlas_info = readtable([Refer_dir 'HCP-MMP1_UniqueRegionList.csv']);



%% read BB gradient and moments, From enigmatools
BB_gradient = load([Refer_dir 'enigma/bb_gradient_glasser_360.csv'])';
BB_moments = load([Refer_dir 'enigma/bb_moments_glasser_360.csv'])';


%% read Economo-Koskinas cytoarchitectonics, From enigmatools
tmp_suffix = '_fsa5'; % _fsa5 _conte69
EKc_parcels = load([Refer_dir 'enigma/economo_koskinas' tmp_suffix '.csv']);
glasser_parcels = load([Refer_dir 'enigma/glasser_360' tmp_suffix '.csv']);
tmpEKc = [EKc_parcels glasser_parcels];
[EKc, EKcw] = subfun_disposalEKcData(tmpEKc);


%% read Yeo2011 function network, From https://github.com/ThomasYeoLab/CBIG/blob/master/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/
Yeo_FN_fsa5 = load([Refer_dir 'Yeo2011_1000subjects_clusters_7_ref.mat']);
Yeo_FN_parcels = [Yeo_FN_fsa5.lh_labels; Yeo_FN_fsa5.rh_labels];
glasser_parcels = load([Refer_dir 'enigma/glasser_360_fsa5.csv']);
tmpYeo = [Yeo_FN_parcels glasser_parcels];
[Yeo_FN, Yeo_FNw] = subfun_disposalYeoFNData(tmpYeo);


%% read Mito data
MitoData = load([Group_dir FileName_Mito]);
zMitoData = zscore(MitoData);


%% read all other data maps, From neuromaps
all_other_maps = readtable([Refer_dir 'all_other_data.csv']);
tmp = table2array(all_other_maps(:,5:end));
tmpMaps = {};
for i = 1:size(tmp,1)
    tmpMaps{i} = tmp(i,:)';
end
tmp = table(tmpMaps', 'VariableNames',{'maps'});
all_other_maps = [all_other_maps(:,1:4) tmp];
all_other_maps([1 5 7 9 11 13 57],:) = [];

all_other_maps_details_des = readtable([Refer_dir 'all_other_data_neuromaps_biological_meaning.xlsx']);
aom_cholinergic_system = [1 3 72 ... % Vesicular acetylcholine transporter (VAChT) – marker of cholinergic terminals
                          29]; % α4β2 nicotinic acetylcholine receptor – cholinergic system
aom_dopaminergic_system = [31 ... % D1 dopamine receptor – dopaminergic system
                           30 62 70 ... % D2/D3 dopamine receptor – dopaminergic system
                           2 37 ... % D2/D3 dopamine receptor antagonist – dopaminergic system
                           13 63]; % Dopamine transporter (DAT) – dopaminergic system
aom_serotonergic_system = [6 68 ... % 5-HT1A receptor – serotonergic system
                           4 16 67 ... % 5-HT1B receptor – serotonergic system
                           5 65 ... % 5-HT2A receptor – serotonergic system
                           52 ... % 5-HT2C receptor – serotonergic system
                           8 ... % 5-HT4 receptor – serotonergic system
                           53 ... % 5-HT6 receptor – serotonergic system
                           7 14 66]; % Serotonin transporter (SERT) – serotonergic system
aom_noradrenergic_system = [10 28]; % Norepinephrine transporter (NET) – noradrenergic system
aom_glutamatergic_system = [49 ... % mGluR1 – glutamatergic system
                            11 34 61 69 ... % Metabotropic glutamate receptor subtype 5 (mGluR5) – glutamatergic system
                            18]; % NMDA receptor – glutamatergic system
aom_GABAergic_system = [12 36 51]; % GABA-A receptor – GABAergic system
aom_histaminergic_system = [17]; % Histamine H3 receptor – histaminergic system
aom_opioidergic_system = [32 73 ... % μ-opioid receptor – opioidergic system
                          74]; % Kappa opioid receptor – opioidergic system



                      
                      
%% read Lipidome maps
Lipid_HRMS = load([Refer_dir '2024NC_LipidAtlas/HRMS_Lipidome_glasserMap.mat']).HRMS_Lipidome;
Lipid_MRM = load([Refer_dir '2024NC_LipidAtlas/MRM_Lipidome_glasserMap.mat']).MRM_Lipidome;







%% ======================= step 01-00 for disposal all data
% input all the basic data

%% set parameters
% input parameters
subjs = {'S01','S02','S03','S04','S05','S06','S07','S08','S09','S10',...
         'S11','S12','S13','S14','S15','S16','S17','S18','S19','S20',...
         'S21','S22','S23','S24','S25','S26','S27','S28','S29','S30'};
which_subjs = [1:3 6 8:9 13 15:17 19:22];
num_layer = 6;
atlas = 'glasser360';
ifloadGO = 1;
signal = 'CBF';
FileName_Mito = ['Mito_' atlas '_bil.1D'];

% basic parameters
Refer_dir = '../../reference/';
Bigbrain_dir = [Refer_dir 'BigBrain/'];
Data_dir = '../Data/';
Group_dir = [Data_dir 'group/'];
CellType_dir = [Refer_dir 'CellType/'];
FileName_BigBrain = ['wholeBB_' atlas '_' num2str(num_layer) 'layers_LaminarProfile.mat'];



%% read scan and BigBrain data
% read scan data s1
scan_suffix = '_s1';
Data_scan = {};
Data_Val = [];
Data_ValLayer = {}; Data_ValLayer{num_layer} = [];
Data_SNR = [];
Data_CNRm = {};
for i = which_subjs
    subj = subjs{i};
    tmp = load([Data_dir subj '/LaminarProfile_' num2str(num_layer) 'points_' atlas '_' signal scan_suffix '.mat']);
    Data_scan = {Data_scan{:}, tmp.LaminarProfile};
    Data_Val = [Data_Val tmp.MRIsignal];
    for k = 1:num_layer
        Data_ValLayer{k} = [Data_ValLayer{k} tmp.LaminarProfile(:,k)];
    end
    Data_SNR = [Data_SNR tmp.SNR];
    Data_CNRm{size(Data_SNR,2)} = tmp.CNRm;
end

Data_Val_bil = Data_Val;
Data_Val_uni = subfun_merge_bil2uni(Data_Val_bil')';
Data_ValLayer_bil = Data_ValLayer;
Data_ValLayer_uni = {};
for i = 1:num_layer
    Data_ValLayer_uni{i} = subfun_merge_bil2uni(Data_ValLayer_bil{i}')';
end

Data_scan_bil = Data_scan';
Data_scan_uni = {};
for i = 1:numel(Data_scan_bil)
    Data_scan_uni{i} = subfun_merge_bil2uni(Data_scan_bil{i}')';
end


Data1_Val_bil = Data_Val_bil;
Data1_Val_uni = Data_Val_uni;
Data1_scan_bil = Data_scan_bil;
Data1_scan_uni = Data_scan_uni;



% read scan data s2
scan_suffix = '_s2';
Data_scan = {};
Data_Val = [];
Data_ValLayer = {}; Data_ValLayer{num_layer} = [];
Data_SNR = [];
Data_CNRm = {};
for i = which_subjs
    subj = subjs{i};
    tmp = load([Data_dir subj '/LaminarProfile_' num2str(num_layer) 'points_' atlas '_' signal scan_suffix '.mat']);
    Data_scan = {Data_scan{:}, tmp.LaminarProfile};
    Data_Val = [Data_Val tmp.MRIsignal];
    for k = 1:num_layer
        Data_ValLayer{k} = [Data_ValLayer{k} tmp.LaminarProfile(:,k)];
    end
    Data_SNR = [Data_SNR tmp.SNR];
    Data_CNRm{size(Data_SNR,2)} = tmp.CNRm;
end

Data_Val_bil = Data_Val;
Data_Val_uni = subfun_merge_bil2uni(Data_Val_bil')';
Data_ValLayer_bil = Data_ValLayer;
Data_ValLayer_uni = {};
for i = 1:num_layer
    Data_ValLayer_uni{i} = subfun_merge_bil2uni(Data_ValLayer_bil{i}')';
end

Data_scan_bil = Data_scan';
Data_scan_uni = {};
for i = 1:numel(Data_scan_bil)
    Data_scan_uni{i} = subfun_merge_bil2uni(Data_scan_bil{i}')';
end


Data2_Val_bil = Data_Val_bil;
Data2_Val_uni = Data_Val_uni;
Data2_scan_bil = Data_scan_bil;
Data2_scan_uni = Data_scan_uni;





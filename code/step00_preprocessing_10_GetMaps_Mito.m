clc
clear all
close all
set(0,'defaultfigurecolor',[1 1 1]);% set figure background color

%% set parameters
% input parameters
atlas = 'glasser360';
file_readname = 'MitoData_HCPMMP1.1D';
MitoNames = {'CI' 'CII' 'CIV' 'MitoD' 'MRC' 'TRC'};

% basic parameters
Data_dir = '/Users/guofanhua/Desktop/gfh/work/StandardBrainTemplateAndAtlas/Mitochondrial_oxidative_phosphorylation/MNI/data/';
Out_dir = '/Users/guofanhua/Desktop/gfh/work/experiment/ASL_Mesoscopic2025/code/Data/group/';

%% read data
% read data
Data = load([Data_dir file_readname]);
Data = Data(:,4:end);


%% disposal
MitoData = [];
for i = 1:360
    tmpData = Data(Data(:,1)==i,2:end);
    tmp = [];
    for k = 1:6
        tmp = [tmp mean(tmpData(tmpData(:,k)~=0,k))];
    end
    MitoData = [MitoData; tmp];
end

MitoData_bil = MitoData;
MitoData_uni = (MitoData(1:end/2,:)+MitoData(end/2+1:end,:))/2;
MitoData_uni = [MitoData_uni; MitoData_uni];

%% save
save([Out_dir 'Mito_' atlas '_uni.1D'],'MitoData_uni','-ascii');
save([Out_dir 'Mito_' atlas '_bil.1D'],'MitoData_bil','-ascii');

name = {'CI' 'CII' 'CIV' 'MitoD' 'MRC' 'TRC'};
for i = 1:6
    tmp = MitoData_uni(:,i);
    save([Out_dir 'Mito_' name{i} '_uni.1D'],'tmp','-ascii');
    
    tmp = MitoData_bil(:,i);
    save([Out_dir 'Mito_' name{i} '_bil.1D'],'tmp','-ascii');
end





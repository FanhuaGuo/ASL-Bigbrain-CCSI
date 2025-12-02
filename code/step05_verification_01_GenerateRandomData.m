clc
clear all
close all
set(0,'defaultfigurecolor',[1 1 1]);% set figure background color
addpath(genpath('/Users/guofanhua/Desktop/gfh/tools/MatlabTools/BrainSpace-latest/matlab'))

%% ======================= step 1 for disposal all data

%% input data
step01_CBFandCCSImaps_00_init;
step01_CBFandCCSImaps_00_ComputeBasicResults;


%% genarate
rng(1234);
CCSI_bil = [];
CCSI_uni = [];
TotalNumber = size(Data_PScs_bil,2);
RepeatTimes = 20;
OutDir = '../Data/group/verification';

% 5 subjects
PickNum = 5;
for i = 1:RepeatTimes
    tmpPick = randperm(TotalNumber);
    CCSI_bil = mean(Data_PScs_bil(:,tmpPick(1:PickNum)),2);
    CCSI_uni = mean(Data_PScs_uni(:,tmpPick(1:PickNum)),2);
    CCSI_uni = [CCSI_uni; CCSI_uni];
    
    suffix = '_bil.1D';
    save([OutDir '/CCSI_' num2str(PickNum) '_' num2str(i) suffix],'CCSI_bil','-ascii');
    
    suffix = '_uni.1D';
    save([OutDir '/CCSI_' num2str(PickNum) '_' num2str(i) suffix],'CCSI_uni','-ascii');
end


% 10 subjects
PickNum = 10;
for i = 1:RepeatTimes
    tmpPick = randperm(TotalNumber);
    CCSI_bil = mean(Data_PScs_bil(:,tmpPick(1:PickNum)),2);
    CCSI_uni = mean(Data_PScs_uni(:,tmpPick(1:PickNum)),2);
    CCSI_uni = [CCSI_uni; CCSI_uni];
    
    suffix = '_bil.1D';
    save([OutDir '/CCSI_' num2str(PickNum) '_' num2str(i) suffix],'CCSI_bil','-ascii');
    
    suffix = '_uni.1D';
    save([OutDir '/CCSI_' num2str(PickNum) '_' num2str(i) suffix],'CCSI_uni','-ascii');
end


% 15 subjects
PickNum = 15;
for i = 1:RepeatTimes
    tmpPick = randperm(TotalNumber);
    CCSI_bil = mean(Data_PScs_bil(:,tmpPick(1:PickNum)),2);
    CCSI_uni = mean(Data_PScs_uni(:,tmpPick(1:PickNum)),2);
    CCSI_uni = [CCSI_uni; CCSI_uni];
    
    suffix = '_bil.1D';
    save([OutDir '/CCSI_' num2str(PickNum) '_' num2str(i) suffix],'CCSI_bil','-ascii');
    
    suffix = '_uni.1D';
    save([OutDir '/CCSI_' num2str(PickNum) '_' num2str(i) suffix],'CCSI_uni','-ascii');
end


% 20 subjects
PickNum = 20;
for i = 1:RepeatTimes
    tmpPick = randperm(TotalNumber);
    CCSI_bil = mean(Data_PScs_bil(:,tmpPick(1:PickNum)),2);
    CCSI_uni = mean(Data_PScs_uni(:,tmpPick(1:PickNum)),2);
    CCSI_uni = [CCSI_uni; CCSI_uni];
    
    suffix = '_bil.1D';
    save([OutDir '/CCSI_' num2str(PickNum) '_' num2str(i) suffix],'CCSI_bil','-ascii');
    
    suffix = '_uni.1D';
    save([OutDir '/CCSI_' num2str(PickNum) '_' num2str(i) suffix],'CCSI_uni','-ascii');
end




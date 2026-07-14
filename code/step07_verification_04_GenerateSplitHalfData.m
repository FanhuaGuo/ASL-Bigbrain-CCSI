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


%% genarate split-half data
rng(1234);
CCSI_bil = [];
CCSI_uni = [];
TotalNumber = size(Data_PScs_bil_correct,2);
RepeatTimes = 20;
PickNum = 15;
OutDir = '../Data/group/SplitHalf';
% mkdir(OutDir)

% run
r = [];
for i = 1:RepeatTimes
    tmpPick = randperm(TotalNumber);
    CCSI_bil1 = mean(Data_PScs_bil_correct(:,tmpPick(1:PickNum)),2);
    CCSI_uni1 = mean(Data_PScs_uni(:,tmpPick(1:PickNum)),2);
    CCSI_uni1 = [CCSI_uni1; CCSI_uni1];
    
    CCSI_bil2 = mean(Data_PScs_bil_correct(:,tmpPick(PickNum+1:end)),2);
    CCSI_uni2 = mean(Data_PScs_uni(:,tmpPick(PickNum+1:end)),2);
    CCSI_uni2 = [CCSI_uni2; CCSI_uni2];
    
    r = [r; corr(CCSI_bil1,CCSI_bil2)];
    
    suffix = '_bil.1D';
    save([OutDir '/CCSI1_' num2str(i) suffix],'CCSI_bil1','-ascii');
    suffix = '_uni.1D';
    save([OutDir '/CCSI1_' num2str(i) suffix],'CCSI_uni1','-ascii');
    suffix = '_bil.1D';
    save([OutDir '/CCSI2_' num2str(i) suffix],'CCSI_bil2','-ascii');
    suffix = '_uni.1D';
    save([OutDir '/CCSI2_' num2str(i) suffix],'CCSI_uni2','-ascii');
end
















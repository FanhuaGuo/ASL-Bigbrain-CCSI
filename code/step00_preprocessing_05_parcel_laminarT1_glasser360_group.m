

%%
clc
clear all
close all
set(0,'defaultfigurecolor',[1 1 1]);% set figure background color

% set parameters
% input parameters

subjs = {'S01' 'S02' 'S03' 'S04' 'S05' 'S06' 'S07' 'S08' 'S09' 'S13' 'S15' 'S16' 'S17' 'S18' 'S19' 'S20' 'S21' 'S22' 'S32' 'S36' 'S37' 'S38' 'S39' 'S40' 'S41' 'S42' 'S43' 'S44' 'S45' 'S46'};
atlasFrom = 2;  % 1:MNI152;  2:std141

for kkkk = 1:numel(subjs)
subj = subjs{kkkk}
% subj = 'S01';
num_layer = 50;

switch atlasFrom
    case 1
        atlas = 'glasser360';
        atlas_readname = ['Data.HCPMMP1.depth.T1.1D'];
    case 2
        atlas = 'glasser360std141';
        atlas_readname = ['Data.std141HCPMMP1.depth.T1.1D'];
end
% basic parameters
Data_dir = '/Users/guofanhua/Desktop/gfh/work/experiment/ASL_Mesoscopic2025/code/Data';


% read scan data and compute laminar profile
% read data
Data = load([Data_dir '/' subj '/' atlas_readname]);
Data(Data(:,5)<=0,:) = [];  % extract GM voxels
Data(Data(:,5)>1,:) = [];
Data = Data(:,4:end);





% compute laminar profile
LaminarProfileT1w = [];
LaminarProfileT1v = [];
LaminarProfile_voxelsNum = [];
for roi = 1:360
    % extract ROI i data
    tmpData = Data(Data(:,1)==roi,2:end);
    
    % compute laminar profile
    tmpLPw = [];
    tmpLPv = [];
    tmpLPVN = [];
    for i = 1:num_layer
        tmpROI = tmpData(:,1)>(i-1)*1/num_layer & tmpData(:,1)<=i*1/num_layer;
        tmpLPw = [tmpLPw mean(tmpData(tmpROI,2),1)];
        tmpLPv = [tmpLPv mean(tmpData(tmpROI,3),1)];
        tmpLPVN = [tmpLPVN sum(tmpROI)];
    end
    LaminarProfileT1w = [LaminarProfileT1w; tmpLPw];
    LaminarProfileT1v = [LaminarProfileT1v; tmpLPv];
    LaminarProfile_voxelsNum = [LaminarProfile_voxelsNum; tmpLPVN];
end

% save
save([Data_dir '/' subj '/LaminarProfile_' num2str(num_layer) 'points_' atlas '_T1wv.mat'],'LaminarProfileT1w','LaminarProfileT1v','LaminarProfile_voxelsNum');





end



clc
clear all
close all
set(0,'defaultfigurecolor',[1 1 1]);% set figure background color

%% set parameters
% input parameters
subj = 'S01';
session = 1; % s1 or s2
num_layer = 20;
atlasFrom = 2;  % 1:MNI152;  2:std141

switch atlasFrom
    case 1
        atlas = 'glasser360';
        atlas_readname = ['Data.HCPMMP1.depth.CBF.s' num2str(session) '.1D'];
    case 2
        atlas = 'glasser360std141';
        atlas_readname = ['Data.std141HCPMMP1.depth.CBF.s' num2str(session) '.1D'];
end
% basic parameters
Data_dir = '/Users/guofanhua/Desktop/gfh/work/experiment/ASL_Mesoscopic2025/code/Data';


%% read scan data and compute laminar profile
% read data
Data = load([Data_dir '/' subj '/' atlas_readname]);
Data(Data(:,5)<=0,:) = [];  % extract GM voxels
Data(Data(:,5)>1,:) = [];
Data = Data(:,4:end);
for i = 3:size(Data,2)    % remove extreme voxels
    Data(Data(:,i)==0,:) = [];
    Data(Data(:,i)>200,:) = [];
    Data(Data(:,i)<-200,:) = [];
end

% compute laminar profile
LaminarProfile = [];
LaminarProfile_voxelsNum = [];
MRIsignal = [];
SNR = [];
total_mu = [];
total_sigma = [];
for roi = 1:360
    % extract ROI i data
    tmpData = Data(Data(:,1)==roi,2:end);
    
    % get avg CBF
    MRIsignal = [MRIsignal; mean(mean(tmpData(:,2:end),2),1)];
    
    % compute SNR in ROI i:  split-half method
    tt = tmpData(:,2:end);
    switch size(tt,2)
        case 2  % σ = | CBF1 - CBF2 | / 2
            tmp = abs(tt(:,1)-tt(:,2))/2;
            sigma = median(tmp);
            mu = mean(tt(:));
            tmpSNR = mu/sigma;
            SNR = [SNR; tmpSNR];
        case 3  % σ = median(|mean(CBF1,CBF2)-CBF3|,|mean(CBF1,CBF3)-CBF2|,|mean(CBF2,CBF3)-CBF1|) / sqrt(4.5)
            tmp = [abs(tt(:,1)/2 + tt(:,2)/2 - tt(:,3)) ...
                   abs(tt(:,1)/2 + tt(:,3)/2 - tt(:,2)) ...
                   abs(tt(:,2)/2 + tt(:,3)/2 - tt(:,1))]/sqrt(4.5);
            sigma = median(median(tmp,2),1);
            mu = mean(tt(:));
            tmpSNR = mu/sigma;
            SNR = [SNR; tmpSNR];
        otherwise  % n>=4: σ = | mu~CBFodd - mu~CBFeven | / sqrt( N~total * ( 1/N~odd + 1/N~even) )
            odd = tt(:,1:2:end);
            even = tt(:,2:2:end);
            tmp = abs(mean(odd,2)-mean(even,2))/sqrt(size(tt,2)*(1/size(odd,2)+1/size(even,2)));
            sigma = median(tmp);
            mu = mean(tt(:));
            tmpSNR = mu/sigma;
            SNR = [SNR; tmpSNR];
    end
    total_mu = [total_mu; mu];
    total_sigma = [total_sigma; sigma];
    
    % compute laminar profile
    tmpLP = [];
    tmpLPVN = [];
    for i = 1:num_layer
        tmpROI = tmpData(:,1)>(i-1)*1/num_layer & tmpData(:,1)<=i*1/num_layer;
        tmpLP = [tmpLP mean(mean(tmpData(tmpROI,2:end),2),1)];
        tmpLPVN = [tmpLPVN sum(tmpROI)];
    end
    LaminarProfile = [LaminarProfile; tmpLP];
    LaminarProfile_voxelsNum = [LaminarProfile_voxelsNum; tmpLPVN];
end


% compute CNR matrix
CNRm = zeros(360);
for i = 1:360
    for j = 1:360
        CNRm(i,j) = (total_mu(i) - total_mu(j)) / sqrt(0.5 * (total_sigma(i)^2 + total_sigma(j)^2));
    end
end


%% draw
figure; hold on;
for i = 1:180
    plot(LaminarProfile(i,:));
end

figure; hold on;
for i = 1:180
    plot(LaminarProfile_voxelsNum(i,:));
end


figure;
plot(LaminarProfile(1,:))

figure;
plot(SNR,MRIsignal,'o');

figure;
imagesc(CNRm);


%% save
save([Data_dir '/' subj '/LaminarProfile_' num2str(num_layer) 'points_' atlas '_CBF_s' num2str(session) '.mat'],'LaminarProfile','MRIsignal','SNR','CNRm');







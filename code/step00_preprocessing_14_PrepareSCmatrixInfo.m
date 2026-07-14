clc
clear all
close all
set(0,'defaultfigurecolor',[1 1 1]);% set figure background color
addpath('/Users/guofanhua/Desktop/gfh/tools/MatlabTools/BrainConnectivityToolbox/2019_03_03_BCT')
addpath('/Users/guofanhua/Desktop/gfh/work/experiment/ASL_Mesoscopic2025/reference/2019PNAS_FCSC/code')
addpath('/Users/guofanhua/Desktop/gfh/tools/MatlabTools/BrainSpace-latest/matlab')
addpath('/Users/guofanhua/Desktop/gfh/tools/MatlabTools/ENIGMA-2.0.0/matlab')

%% note
% Fanhua adjusted from Vázquez-Rodríguez B, Suárez L E, Markello R D, et al. Gradients of structure–function tethering across neocortex[J]. Proceedings of the National Academy of Sciences, 2019, 116(42): 21219-21227.
%% set parameters
% input parameters


% basic parameters
Refer_dir = '../../reference/';
Data_dir = '../Data/';
Group_dir = [Data_dir 'group/'];


%% read data
scm = load([Refer_dir 'enigma/strucMatrix_ctx_glasser_360.csv']);  % data from enigmatools

fcm = load([Refer_dir 'enigma/funcMatrix_ctx_glasser_360.csv']);  % data from enigmatools
gm = GradientMaps('kernel','cs','approach','dm','n_components',10,'alignment','procrustes'); % none procrustes
gm = gm.fit(fcm);
fcm_g = gm.gradients{1};  % did not use, I directly used the fcgradient1-10 in neuromaps 

tmp_suffix = '_fsa5'; % _fsa5 _conte69
glasser_parcels = load([Refer_dir 'enigma/glasser_360' tmp_suffix '.csv']);

load('tmpSNRmap.mat');
NTSMTmaps = load('tmpNeurotransmitter.mat').Data1;



%% compute the centroids of all rois in glasser360 atlas and euclidean distance matrix
V = spm_vol('/Users/guofanhua/Desktop/gfh/work/StandardBrainTemplateAndAtlas/MNI_Glasser_HCP_2019_v1.0/HCP_atlas.nii');
Y = spm_read_vols(V);

% 取所有 ROI 标签，去掉背景 0
labels = unique(round(Y(:)));
labels(labels == 0) = [];

centroids = zeros(numel(labels), 3);

for n = 1:numel(labels)
    idx = find(round(Y) == labels(n));
    [i, j, k] = ind2sub(size(Y), idx);

    % 体素中心的平均位置（1-based voxel index）
    vox_center = [mean(i), mean(j), mean(k), 1];

    % 转成世界坐标(mm)
    mm_center = (V.mat * vox_center')';
    centroids(n, :) = mm_center(1:3);
end

euclidean_distance_m = squareform(pdist(centroids));  % euclidean distance
eu = euclidean_distance_m;

%% compute path length (brain connectivity toolbox) and communicability matrix
sc = scm;
sc(sc > 0) = 1;
sc(sc < 0) = 0;

path_length_m = distance_bin(sc);         % path length (brain connectivity toolbox)
communicability_m = fcn_communicability(sc);  % communicability

sp = path_length_m;
co = communicability_m;


%% test regression
clc
mask = CBFsnr_mean_bil>0;
% mask = [1:360]' <= 180;


fc = fcm;
fc = fc(mask,mask);
n = length(fc);

sp = path_length_m(mask,mask);         % path length (brain connectivity toolbox)
co = communicability_m(mask,mask);  % communicability
eu = euclidean_distance_m(mask,mask);  % euclidean distance

% initialize node-wise R-square values for resolution ii
rsq = zeros(n, 1);
for jj = 1:n
    % define fc response (y) and sc predictors (x_i)
    y = fc(:, jj);

    x1 = sp(:, jj);
    x2 = co(:, jj);
    x3 = eu(:, jj);

    % standardize predictors
    x = zscore([x1, x2, x3]);

    % fit multiple regression (OLS, main effects only)
    % exclude self-connections for all variables
    % N.B., `fitlm` adds an intercept by default
    lm = fitlm(x, y, 'Exclude', jj);
    % record adjusted R-square for parcellation ii, node jj
    rsq(jj) = lm.Rsquared.Adjusted;

    % let me know how we're doing :)
%     fprintf('node %i out of %i done\n', jj, n)
end
tt = zeros(360,1);
tt(mask) = rsq;
rsq = tt;
corr(-fcm_g(mask,1),rsq(mask))



DataName = ['CBF_Mean_PScs_L12'];
CCSI = load([Group_dir DataName  '_bil.1D']);
CCSI = CCSI(mask);
% CCSI = abs(CCSI - CCSI');
% CCSI = abs(CCSI * CCSI');
% CCSI = CCSI(mask,mask);
BB_gradient = load([Refer_dir 'enigma/bb_gradient_glasser_360.csv'])';
BB_gradient = BB_gradient(mask);
NTSm = corr(zscore(NTSMTmaps)');
NTSm = NTSm(mask,mask);

rsq_ccsi = zeros(n, 1);
for jj = 1:n
    % define fc response (y) and sc predictors (x_i)
    y = fc(:, jj);

    x1 = sp(:, jj);
    x2 = co(:, jj);
    x3 = eu(:, jj);
%     x4 = CCSI;
%     x4 = CCSI(:,jj);
%     x4 = BB_gradient;
    x4 = NTSm(:, jj);

    % standardize predictors
    x = zscore([x1, x2, x3, x4]);

    % fit multiple regression (OLS, main effects only)
    % exclude self-connections for all variables
    % N.B., `fitlm` adds an intercept by default
    lm = fitlm(x, y, 'Exclude', jj);
    % record adjusted R-square for parcellation ii, node jj
    rsq_ccsi(jj) = lm.Rsquared.Adjusted;  % Adjusted Ordinary

    % let me know how we're doing :)
%     fprintf('node %i out of %i done\n', jj, n)
end
tt = zeros(360,1);
tt(mask) = rsq_ccsi;
rsq_ccsi = tt;

figure; hold on;
plot(-fcm_g(mask,1),rsq(mask),'bo')
plot(-fcm_g(mask,1),rsq_ccsi(mask),'ro')

figure;
plot(-fcm_g(mask,1),rsq_ccsi(mask)-rsq(mask),'bo')


figure; hold on;
plot(zscore(rsq_ccsi(mask)-rsq(mask)),'b-');
plot(zscore(-fcm_g(mask,1)),'r-');


corr(-fcm_g(mask,1),rsq_ccsi(mask))

corr(-fcm_g(mask,1),rsq_ccsi(mask)-rsq(mask))


%% plotting
tmpOut = zeros(size(glasser_parcels));
for i = 1:360
    tmpOut(glasser_parcels==i) = rsq_ccsi(i)-rsq(i);
end
figure;
plot_cortical(tmpOut, ...
    'surface_name', 'fsa5', ...
    'label_text', 'rsq ccsi - rsq');   % colormap我们自己定义，不用默认


tmpOut = zeros(size(glasser_parcels));
for i = 1:360
    tmpOut(glasser_parcels==i) = rsq(i);
end
figure;
plot_cortical(tmpOut, ...
    'surface_name', 'fsa5', ...
    'label_text', 'rsq');   % colormap我们自己定义，不用默认


tmpOut = zeros(size(glasser_parcels));
for i = 1:360
    tmpOut(glasser_parcels==i) = rsq_ccsi(i);
end
figure;
plot_cortical(tmpOut, ...
    'surface_name', 'fsa5', ...
    'label_text', 'rsq ccsi');   % colormap我们自己定义，不用默认


tmpOut = zeros(size(glasser_parcels));
for i = 1:360
    tmpOut(glasser_parcels==i) = CCSI(i);
end
figure;
plot_cortical(tmpOut, ...
    'surface_name', 'fsa5', ...
    'label_text', 'ccsi');   % colormap我们自己定义，不用默认



%% test
close all
jj = 30;
% define fc response (y) and sc predictors (x_i)
y = fc(:, jj);

x1 = sp(:, jj);
x2 = co(:, jj);
x3 = eu(:, jj);
x4 = CCSI;

% standardize predictors
x = zscore([x1, x2, x3, x4]);

% fit multiple regression (OLS, main effects only)
% exclude self-connections for all variables
% N.B., `fitlm` adds an intercept by default
lm = fitlm(x, y, 'Exclude', jj);
% record adjusted R-square for parcellation ii, node jj
lm.Rsquared.Ordinary

% let me know how we're doing :)
%     fprintf('node %i out of %i done\n', jj, n)

figure; hold on;
plot(x(:,1),'k-');
plot(x(:,2),'b-');
plot(x(:,3),'g-');
plot(x(:,4));
plot(zscore(y),'r-');





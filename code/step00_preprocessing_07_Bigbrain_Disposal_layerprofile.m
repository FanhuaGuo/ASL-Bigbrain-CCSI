clc
clear all
close all
set(0,'defaultfigurecolor',[1 1 1]);% set figure background color

% compare laminar profile of perfusion signals and cell body density
%% set parameters
% basic parameters
NumROIs = 360;
Numlayers = 6;
AnalysisDir = '/Users/guofanhua/Desktop/gfh/work/StandardBrainTemplateAndAtlas/BigBrain_LaminarProfile/whole_SUMA/LaminarProfile_HCP_MMP1p0/';
SaveName = '/Users/guofanhua/Desktop/gfh/work/StandardBrainTemplateAndAtlas/BigBrain_LaminarProfile/Results/wholeBB_6layers_LaminarProfile.mat';
ifsave = 0;

%% read and disposal
LP = [];
for i = 1:NumROIs
    tmpLP = [];
    for j = 1:Numlayers
        tmp = load([AnalysisDir 'ROI' num2str(i) '_depth' num2str(j) '.1D']);
        tmpLP = [tmpLP tmp];
    end
    LP = [LP; tmpLP];
end

LP = 65535-LP;
BB_LaminarProfile = LP;

%% draw

figure; hold on;
for i = 1:NumROIs
    plot(LP(i,:));
end

%% save
if ifsave
    save(SaveName,'BB_LaminarProfile');
end


%% draw V1
Numlayers = 50;
AnalysisDir = '../reference/BigBrain/LaminarProfile_HCP_MMP1p0_V1/';
V1LP = [];
for j = 3:Numlayers-2
    tmp = load([AnalysisDir 'V1_depth' num2str(j) '.1D']);
    V1LP = [V1LP; tmp];
end
V1LP = 65535-V1LP;

figure; hold on;
plot(3:Numlayers-2,V1LP,'k-','linewidth',4);
set(gca,'xlim',[0 Numlayers+1],'xTick',[1 Numlayers],'XTickLabels',{'WM','CSF'},'FontSize',16,'linewidth',2);
set(gca,'yTick',[],'FontSize',16,'linewidth',2);
ylabel('cell-body staining intensity','FontSize',16,'linewidth',10);
box off





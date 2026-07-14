clc
clear all
close all
set(0,'defaultfigurecolor',[1 1 1]);% set figure background color

% compare laminar profile of perfusion signals and cell body density
%% set parameters
% basic parameters
NumROIs = 360;
Numlayers = 15;
AnalysisDir = ['/Users/guofanhua/Desktop/gfh/work/StandardBrainTemplateAndAtlas/BigBrain_LaminarProfile/whole_SUMA/LaminarProfile_HCP_MMP1_layer' num2str(Numlayers) '/'];
SavePath = '/Users/guofanhua/Desktop/gfh/work/StandardBrainTemplateAndAtlas/BigBrain_LaminarProfile/Results/';
% SaveName = ['wholeBB_' num2str(Numlayers) 'layers_LaminarProfile.mat'];
ifsave = 1;

%% read and disposal
LP = [];
for i = 1:NumROIs
    tmpLP = [];
    for j = 1:Numlayers
        tmp = load([AnalysisDir 'ROI' num2str(i) '_depth' num2str(j) '.1D']);
        if isempty(tmp) tmp = 0; end
        tmpLP = [tmpLP tmp];
    end
    LP = [LP; tmpLP];
end

LP = 65535-LP;
BB_LaminarProfile = LP;


%% correct
Label_NoSignal = [];
for i = 1:size(LP,2)
    if mean(LP(:,i)) > 60000
        Label_NoSignal = [Label_NoSignal; 0];
    else
        Label_NoSignal = [Label_NoSignal; 1];
    end
end
outNumlayers = sum(Label_NoSignal);
BB_LaminarProfile = BB_LaminarProfile(:,Label_NoSignal>0);

%% draw
figure; hold on;
for i = 1:NumROIs
    plot(BB_LaminarProfile(i,:));
end

%% save
if ifsave
    SaveName = ['wholeBB_glasser360_' num2str(outNumlayers) 'layers_LaminarProfile.mat'];
    save([SavePath SaveName],'BB_LaminarProfile');
end


%% draw V1
drawV1 = 1;
if drawV1
    Numlayers = 50;
    AnalysisDir = '/Users/guofanhua/Desktop/gfh/work/StandardBrainTemplateAndAtlas/BigBrain_LaminarProfile/whole_SUMA/LaminarProfile_HCP_MMP1p0_V1/';
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
end





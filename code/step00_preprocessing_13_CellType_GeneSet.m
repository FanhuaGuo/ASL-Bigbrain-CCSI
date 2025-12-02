clc
clear all
close all
set(0,'defaultfigurecolor',[1 1 1]);% set figure background color


%% ===================================================
% ===============2025 PB perfusion====================
% ====================================================
ifrun = 0;
if ifrun
%% set parameters
which_gene_data = 2;  % 1:mine AHBA;  2:NN AHBA DME

switch which_gene_data
    case 1
        % input parameters
        FileName = '2025PlosB_perfusion.xlsx';
        atlas = 'glasser360';
        % basic parameters
        Data_dir = '/Users/guofanhua/Desktop/gfh/work/experiment/ASL_Mesoscopic2025/reference/CellType/';
        Allen_dir = '/Users/guofanhua/Desktop/gfh/work/StandardBrainTemplateAndAtlas/AllenBrain/microarray/';
        FileName_GeneExpression = ['GeneExpression_' atlas '.csv'];
        SaveName = [Data_dir 'CellType_2025PlosB_' atlas '.mat'];
    case 2
        % input parameters
        FileName = '2025PlosB_perfusion.xlsx';
        atlas = 'glasser360';
        % basic parameters
        Data_dir = '/Users/guofanhua/Desktop/gfh/work/experiment/ASL_Mesoscopic2025/reference/CellType/';
        Allen_dir = '/Users/guofanhua/Desktop/gfh/work/experiment/ASL_Mesoscopic2025/reference/2024NNcode_AHBA_gradients-master/outputs/';
        FileName_GeneExpression = 'expression/hcp_3d_gxrm.csv';
        SaveName = [Data_dir 'CellType_2025PlosB_' atlas '_NN.mat'];
end

%% read gene data
Data_AHBA_GE = readtable([Allen_dir FileName_GeneExpression]);

tmpGeneID = Data_AHBA_GE.Var1(2:end);
for i = 1:numel(tmpGeneID)
    if tmpGeneID{i}(1) >= '0' && tmpGeneID{i}(1) <= '9'
        tmpGeneID{i} = subfun_replace_dateGenes(tmpGeneID{i});
    end
end
GeneID = tmpGeneID;

%% read cell type data and extract
Data_CellType_Orig = readtable([Data_dir FileName]);

clear CT_Genes ASCII_Genes Coord_ASCIIGenes
CT_Genes = {};
ASCII_Genes = {};   
for i = 1:26 ASCII_Genes{i} = {' '}; end
Coord_ASCIIGenes = {};   Coord_ASCIIGenes{26} = [];
for i = 1:size(Data_CellType_Orig,1)
    tmp = Data_CellType_Orig.cGenes{i};
    tmps = regexp(tmp,' ');
    tmpc = regexp(tmp,',');
    tt = tmp(1:tmpc(1)-1);
    ttt = ASCII_Genes{double(tt(1))-double('A')+1};
    ASCII_Genes{double(tt(1))-double('A')+1} = {ttt{:}, tt};
    Coord_ASCIIGenes{double(tt(1))-double('A')+1} = [Coord_ASCIIGenes{double(tt(1))-double('A')+1}; i 1];
    tmpGenes = {tt};
    for k = 2:numel(tmpc)
        tt = tmp(tmps(k-1)+1:tmpc(k)-1);
        ttt = ASCII_Genes{double(tt(1))-double('A')+1};
        ASCII_Genes{double(tt(1))-double('A')+1} = {ttt{:}, tt};
        Coord_ASCIIGenes{double(tt(1))-double('A')+1} = [Coord_ASCIIGenes{double(tt(1))-double('A')+1}; i k];
        tmpGenes = {tmpGenes{:}, tt};
    end
    tt = tmp(tmps(end)+1:end);
    ttt = ASCII_Genes{double(tt(1))-double('A')+1};
    ASCII_Genes{double(tt(1))-double('A')+1} = {ttt{:}, tt};
    Coord_ASCIIGenes{double(tt(1))-double('A')+1} = [Coord_ASCIIGenes{double(tt(1))-double('A')+1}; i k];
    tmpGenes = {tmpGenes{:}, tt};
    CT_Genes = {CT_Genes{:},tmpGenes};
end
for i = 1:26 ASCII_Genes{i}(1) = []; end

%% disposal
tmp1 = zeros(size(Data_CellType_Orig,1),1);
tmp2 = {}; 
tmp2{size(Data_CellType_Orig,1)} = [];
Data_CellType = table(Data_CellType_Orig.cLabel, Data_CellType_Orig.cSize, CT_Genes', tmp1, tmp2',...
    'VariableNames',{'Cell', 'OrigSize', 'OrigGenes', 'Size', 'gene_label'});

for i = 1:size(Data_CellType,1) Genes_exist{i} = {' '}; end
for g = 1:numel(GeneID)
    tmp_asc = double(GeneID{g}(1))-double('A')+1;
    for i = 1:numel(ASCII_Genes{tmp_asc})
        switch ASCII_Genes{tmp_asc}{i}
            case GeneID{g}
                tmpCell = Coord_ASCIIGenes{tmp_asc}(i,1);
                Data_CellType.Size(tmpCell) = Data_CellType.Size(tmpCell)+1;
                Data_CellType.gene_label{tmpCell} = [Data_CellType.gene_label{tmpCell}; g];
                Genes_exist{tmpCell} = {Genes_exist{tmpCell}{:}, GeneID{g}};
        end
    end
end
for i = 1:size(Data_CellType,1) Genes_exist{i}(1) = []; end
tmp = table(Genes_exist','VariableNames',{'Genes'});
Data_CellType = [Data_CellType tmp];


%% save
save(SaveName,'Data_CellType');


end



%% ===================================================
% ===============2020 NC & 2024 NN====================
% ====================================================
ifrun = 0;
if ifrun
%% set parameters
which_gene_data = 2;  % 1:mine AHBA;  2:NN AHBA DME

switch which_gene_data
    case 1
        % input parameters
        FileName = '2024NN_2020NC_cellType.csv';
        atlas = 'glasser360';
        % basic parameters
        Data_dir = '/Users/guofanhua/Desktop/gfh/work/experiment/ASL_Mesoscopic2025/reference/CellType/';
        Allen_dir = '/Users/guofanhua/Desktop/gfh/work/StandardBrainTemplateAndAtlas/AllenBrain/microarray/';
        FileName_GeneExpression = ['GeneExpression_' atlas '.csv'];
        SaveName = [Data_dir 'CellType_2020NC_' atlas '_Mine.mat'];
    case 2
        % input parameters
        FileName = '2024NN_2020NC_cellType.csv';
        atlas = 'glasser360';
        % basic parameters
        Data_dir = '/Users/guofanhua/Desktop/gfh/work/experiment/ASL_Mesoscopic2025/reference/CellType/';
        Allen_dir = '/Users/guofanhua/Desktop/gfh/work/experiment/ASL_Mesoscopic2025/reference/2024NNcode_AHBA_gradients-master/outputs/';
        FileName_GeneExpression = 'expression/hcp_3d_gxrm.csv';
        SaveName = [Data_dir 'CellType_2020NC_' atlas '_NN.mat'];
end


%% read gene data
Data_AHBA_GE = readtable([Allen_dir FileName_GeneExpression]);

tmpGeneID = Data_AHBA_GE.Var1(2:end);
for i = 1:numel(tmpGeneID)
    if tmpGeneID{i}(1) >= '0' && tmpGeneID{i}(1) <= '9'
        tmpGeneID{i} = subfun_replace_dateGenes(tmpGeneID{i});
    end
end
GeneID = tmpGeneID;


%% read cell type data and extract
Data_CellType_Orig = readtable([Data_dir FileName]);

clear CT_Genes ASCII_Genes Coord_ASCIIGenes
CT_Genes = {};
ASCII_Genes = {};   
for i = 1:26 ASCII_Genes{i} = {' '}; end
Coord_ASCIIGenes = {};   Coord_ASCIIGenes{26} = [];
for i = 1:size(Data_CellType_Orig,1)
    tmp = Data_CellType_Orig.Genes{i};
    tmpc = regexp(tmp,',');
    tt = tmp(1:tmpc(1)-1);
    if tt(1) >= '0' && tt(1) <= '9'
        tt = subfun_replace_dateGenes(tt);
    end
    ttt = ASCII_Genes{double(tt(1))-double('A')+1};
    ASCII_Genes{double(tt(1))-double('A')+1} = {ttt{:}, tt};
    Coord_ASCIIGenes{double(tt(1))-double('A')+1} = [Coord_ASCIIGenes{double(tt(1))-double('A')+1}; i 1];
    tmpGenes = {tt};
    for k = 2:numel(tmpc)
        tt = tmp(tmpc(k-1)+1:tmpc(k)-1);
        if tt(1) >= '0' && tt(1) <= '9'
            tt = subfun_replace_dateGenes(tt);
        end
        ttt = ASCII_Genes{double(tt(1))-double('A')+1};
        ASCII_Genes{double(tt(1))-double('A')+1} = {ttt{:}, tt};
        Coord_ASCIIGenes{double(tt(1))-double('A')+1} = [Coord_ASCIIGenes{double(tt(1))-double('A')+1}; i k];
        tmpGenes = {tmpGenes{:}, tt};
    end
    tt = tmp(tmpc(end)+1:end);
    if tt(1) >= '0' && tt(1) <= '9'
        tt = subfun_replace_dateGenes(tt);
    end
    ttt = ASCII_Genes{double(tt(1))-double('A')+1};
    ASCII_Genes{double(tt(1))-double('A')+1} = {ttt{:}, tt};
    Coord_ASCIIGenes{double(tt(1))-double('A')+1} = [Coord_ASCIIGenes{double(tt(1))-double('A')+1}; i k];
    tmpGenes = {tmpGenes{:}, tt};
    CT_Genes = {CT_Genes{:},tmpGenes};
end
for i = 1:26 ASCII_Genes{i}(1) = []; end

%% disposal
tmp1 = zeros(size(Data_CellType_Orig,1),1);
tmp2 = {}; 
tmp2{size(Data_CellType_Orig,1)} = [];
Data_CellType = table(Data_CellType_Orig.Type, Data_CellType_Orig.Class, CT_Genes', tmp1, tmp2',...
    'VariableNames',{'Cell', 'Class', 'OrigGenes', 'Size', 'gene_label'});

for i = 1:size(Data_CellType,1) Genes_exist{i} = {' '}; end
for g = 1:numel(GeneID)
    tmp_asc = double(GeneID{g}(1))-double('A')+1;
    for i = 1:numel(ASCII_Genes{tmp_asc})
        switch ASCII_Genes{tmp_asc}{i}
            case GeneID{g}
                tmpCell = Coord_ASCIIGenes{tmp_asc}(i,1);
                Data_CellType.Size(tmpCell) = Data_CellType.Size(tmpCell)+1;
                Data_CellType.gene_label{tmpCell} = [Data_CellType.gene_label{tmpCell}; g];
                Genes_exist{tmpCell} = {Genes_exist{tmpCell}{:}, GeneID{g}};
        end
    end
end
for i = 1:size(Data_CellType,1) Genes_exist{i}(1) = []; end
tmp = table(Genes_exist','VariableNames',{'Genes'});
Data_CellType = [Data_CellType tmp];


%% save
save(SaveName,'Data_CellType');


end






%% ===================================================
% ===============2020 Nature, Vascular================
% ====================================================
ifrun = 0;
if ifrun
%% set parameters
which_gene_data = 2;  % 1:mine AHBA;  2:NN AHBA DME

FileName_End = 'Endo.zonation.genes_cluster.txt';
FileName_Per = 'Per.zonation.genes_cluster.txt';
FileName_SMC = 'SMC.zonation.genes_cluster.txt';
atlas = 'glasser360';
Input_dir = '/Users/guofanhua/Desktop/gfh/work/experiment/ASL_Mesoscopic2025/reference/2022Nature_VascularAtlas/data/';
Data_dir = '/Users/guofanhua/Desktop/gfh/work/experiment/ASL_Mesoscopic2025/reference/CellType/';
Cell_Name = {'Endo-V1','Endo-V2','Endo-capV1','Endo-capV2','Endo-capA2','Endo-capA1','Endo-A','Endo-VA',...
             'Per2a','Per1v-3','Per1v-2','Per1v-1',...
             'aSMC','vSMC-4','vSMC-3','vSMC-2','vSMC-1'};

switch which_gene_data
    case 1
        % basic parameters
        Allen_dir = '/Users/guofanhua/Desktop/gfh/work/StandardBrainTemplateAndAtlas/AllenBrain/microarray/';
        FileName_GeneExpression = ['GeneExpression_' atlas '.csv'];
        SaveName = [Data_dir 'CellType_2022Nature_' atlas '_Mine.mat'];
    case 2
        % basic parameters
        Allen_dir = '/Users/guofanhua/Desktop/gfh/work/experiment/ASL_Mesoscopic2025/reference/2024NNcode_AHBA_gradients-master/outputs/';
        FileName_GeneExpression = 'expression/hcp_3d_gxrm.csv';
        SaveName = [Data_dir 'CellType_2022Nature_' atlas '_NN.mat'];
end


%% read gene data
Data_AHBA_GE = readtable([Allen_dir FileName_GeneExpression]);

tmpGeneID = Data_AHBA_GE.Var1(2:end);
for i = 1:numel(tmpGeneID)
    if tmpGeneID{i}(1) >= '0' && tmpGeneID{i}(1) <= '9'
        tmpGeneID{i} = subfun_replace_dateGenes(tmpGeneID{i});
    end
end
GeneID = tmpGeneID;

%% read cell type data
Data_End = readtable([Input_dir FileName_End]);
Data_Per = readtable([Input_dir FileName_Per]);
Data_SMC = readtable([Input_dir FileName_SMC]);

tmp = {};
tmp{numel(Cell_Name)} = {};
Data_CellType_Orig = table(Cell_Name',tmp','VariableNames',{'Cell', 'Genes'});

Data = Data_End;
qianzhi = 0;
k = 1;
for i = 1:8
    tmp = {Data.gene{k}};
    while Data.cluster(k) == i
        k = k+1;
        if k > size(Data,1)
            break;
        end
        if Data.cluster(k) == i
            tmp = {tmp{:} Data.gene{k}};
        end
    end
    Data_CellType_Orig.Genes{i+qianzhi} = tmp';
end


Data = Data_Per;
qianzhi = 8;
k = 1;
for i = 1:4
    tmp = {Data.gene{k}};
    while Data.cluster(k) == i
        k = k+1;
        if k > size(Data,1)
            break;
        end
        if Data.cluster(k) == i
            tmp = {tmp{:} Data.gene{k}};
        end
    end
    Data_CellType_Orig.Genes{i+qianzhi} = tmp';
end


Data = Data_SMC;
qianzhi = 12;
k = 1;
for i = 1:5
    tmp = {Data.gene{k}};
    while Data.cluster(k) == i
        k = k+1;
        if k > size(Data,1)
            break;
        end
        if Data.cluster(k) == i
            tmp = {tmp{:} Data.gene{k}};
        end
    end
    Data_CellType_Orig.Genes{i+qianzhi} = tmp';
end

%% extract
clear CT_Genes ASCII_Genes Coord_ASCIIGenes
ASCII_Genes = {};   
for i = 1:26 ASCII_Genes{i} = {' '}; end
Coord_ASCIIGenes = {};   Coord_ASCIIGenes{26} = [];
for i = 1:size(Data_CellType_Orig,1)
    tmp = Data_CellType_Orig.Genes{i};
    tmpGenes = {' '};
    for k = 1:numel(tmp)
        tt = tmp{k};
        if tt(1) >= '0' && tt(1) <= '9'
            tt = subfun_replace_dateGenes(tt);
        end
        ttt = ASCII_Genes{double(tt(1))-double('A')+1};
        ASCII_Genes{double(tt(1))-double('A')+1} = {ttt{:}, tt};
        Coord_ASCIIGenes{double(tt(1))-double('A')+1} = [Coord_ASCIIGenes{double(tt(1))-double('A')+1}; i k];
        tmpGenes = {tmpGenes{:}, tt};
    end
    tmpGenes(1) = [];
end
for i = 1:26 ASCII_Genes{i}(1) = []; end

%% disposal
tmp1 = zeros(size(Data_CellType_Orig,1),1);
tmp2 = {}; 
tmp2{size(Data_CellType_Orig,1)} = [];
Data_CellType = table(Data_CellType_Orig.Cell, Data_CellType_Orig.Genes, tmp1, tmp2',...
    'VariableNames',{'Cell', 'OrigGenes', 'Size', 'gene_label'});

for i = 1:size(Data_CellType,1) Genes_exist{i} = {' '}; end
for g = 1:numel(GeneID)
    tmp_asc = double(GeneID{g}(1))-double('A')+1;
    for i = 1:numel(ASCII_Genes{tmp_asc})
        switch ASCII_Genes{tmp_asc}{i}
            case GeneID{g}
                tmpCell = Coord_ASCIIGenes{tmp_asc}(i,1);
                Data_CellType.Size(tmpCell) = Data_CellType.Size(tmpCell)+1;
                Data_CellType.gene_label{tmpCell} = [Data_CellType.gene_label{tmpCell}; g];
                Genes_exist{tmpCell} = {Genes_exist{tmpCell}{:}, GeneID{g}};
        end
    end
end
for i = 1:size(Data_CellType,1) Genes_exist{i}(1) = []; end
tmp = table(Genes_exist','VariableNames',{'Genes'});
Data_CellType = [Data_CellType tmp];


%% save
save(SaveName,'Data_CellType');


end












%% ===================================================
% ===============2019 Nature, oligo===================
% ====================================================
ifrun = 0;
if ifrun
%% set parameters
which_gene_data = 2;  % 1:mine AHBA;  2:NN AHBA DME

atlas = 'glasser360';
Input_dir = '/Users/guofanhua/Desktop/gfh/work/experiment/ASL_Mesoscopic2025/reference/2019Nature_OligoRNAseq/';
Data_dir = '/Users/guofanhua/Desktop/gfh/work/experiment/ASL_Mesoscopic2025/reference/CellType/';
Cell_Name = {'OPC','COPs','Oligo1','Oligo2','Oligo3',...
             'Oligo4','Oligo5','Oligo6','ImOlGs'};

switch which_gene_data
    case 1
        % basic parameters
        Allen_dir = '/Users/guofanhua/Desktop/gfh/work/StandardBrainTemplateAndAtlas/AllenBrain/microarray/';
        FileName_GeneExpression = ['GeneExpression_' atlas '.csv'];
        SaveName = [Data_dir 'CellType_2019Nature_' atlas '_Mine.mat'];
    case 2
        % basic parameters
        Allen_dir = '/Users/guofanhua/Desktop/gfh/work/experiment/ASL_Mesoscopic2025/reference/2024NNcode_AHBA_gradients-master/outputs/';
        FileName_GeneExpression = 'expression/hcp_3d_gxrm.csv';
        SaveName = [Data_dir 'CellType_2019Nature_' atlas '_NN.mat'];
end


%% read gene data
Data_AHBA_GE = readtable([Allen_dir FileName_GeneExpression]);

tmpGeneID = Data_AHBA_GE.Var1(2:end);
for i = 1:numel(tmpGeneID)
    if tmpGeneID{i}(1) >= '0' && tmpGeneID{i}(1) <= '9'
        tmpGeneID{i} = subfun_replace_dateGenes(tmpGeneID{i});
    end
end
GeneID = tmpGeneID;


%% read cell type data
tmp = {};
tmp{numel(Cell_Name)} = {};
Data_CellType_Orig = table(Cell_Name',tmp','VariableNames',{'Cell', 'Genes'});
for i = 1:numel(Cell_Name)
    tmp = readtable([Input_dir Cell_Name{i} '.csv']);
    Data_CellType_Orig.Genes{i} = tmp.gene;
end


%% extract
clear CT_Genes ASCII_Genes Coord_ASCIIGenes
ASCII_Genes = {};   
for i = 1:26 ASCII_Genes{i} = {' '}; end
Coord_ASCIIGenes = {};   Coord_ASCIIGenes{26} = [];
for i = 1:size(Data_CellType_Orig,1)
    tmp = Data_CellType_Orig.Genes{i};
    tmpGenes = {' '};
    for k = 1:numel(tmp)
        tt = tmp{k};
        if tt(1) >= '0' && tt(1) <= '9'
            tt = subfun_replace_dateGenes(tt);
        end
        ttt = ASCII_Genes{double(tt(1))-double('A')+1};
        ASCII_Genes{double(tt(1))-double('A')+1} = {ttt{:}, tt};
        Coord_ASCIIGenes{double(tt(1))-double('A')+1} = [Coord_ASCIIGenes{double(tt(1))-double('A')+1}; i k];
        tmpGenes = {tmpGenes{:}, tt};
    end
    tmpGenes(1) = [];
end
for i = 1:26 ASCII_Genes{i}(1) = []; end

%% disposal
tmp1 = zeros(size(Data_CellType_Orig,1),1);
tmp2 = {}; 
tmp2{size(Data_CellType_Orig,1)} = [];
Data_CellType = table(Data_CellType_Orig.Cell, Data_CellType_Orig.Genes, tmp1, tmp2',...
    'VariableNames',{'Cell', 'OrigGenes', 'Size', 'gene_label'});

for i = 1:size(Data_CellType,1) Genes_exist{i} = {' '}; end
for g = 1:numel(GeneID)
    tmp_asc = double(GeneID{g}(1))-double('A')+1;
    for i = 1:numel(ASCII_Genes{tmp_asc})
        switch ASCII_Genes{tmp_asc}{i}
            case GeneID{g}
                tmpCell = Coord_ASCIIGenes{tmp_asc}(i,1);
                Data_CellType.Size(tmpCell) = Data_CellType.Size(tmpCell)+1;
                Data_CellType.gene_label{tmpCell} = [Data_CellType.gene_label{tmpCell}; g];
                Genes_exist{tmpCell} = {Genes_exist{tmpCell}{:}, GeneID{g}};
        end
    end
end
for i = 1:size(Data_CellType,1) Genes_exist{i}(1) = []; end
tmp = table(Genes_exist','VariableNames',{'Genes'});
Data_CellType = [Data_CellType tmp];


%% save
save(SaveName,'Data_CellType');


end






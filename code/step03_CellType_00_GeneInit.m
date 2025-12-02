
%% ======================= step 01-00 for disposal all data
% input all the basic data

%% set parameters
% input parameters
FileName_CellType1 = ['CellType_2025PlosB_' atlas '_NN.mat'];
FileName_CellType2 = ['CellType_2022Nature_' atlas '_NN.mat'];
FileName_CellType3 = ['CellType_2019Nature_' atlas '_NN.mat'];

% basic parameters
Allen_dir = [Refer_dir '2024NNcode_AHBA_gradients-master/outputs/'];
FileName_GeneExpression = ['expression/hcp_3d_gxrm.csv'];
CellType_dir = [Refer_dir 'CellType/'];


%% read and disposal Gene data
% read gene data
tmp_GeneExpression = readtable([Allen_dir FileName_GeneExpression]);

tmpGeneID = tmp_GeneExpression.Var1(2:end);
for i = 1:numel(tmpGeneID)
    if tmpGeneID{i}(1) >= '0' && tmpGeneID{i}(1) <= '9'
        tmpGeneID{i} = subfun_replace_dateGenes(tmpGeneID{i});
    end
end
GeneID = tmpGeneID;
GeneExpression_Orig = table2array(tmp_GeneExpression(2:end,2:end));

% read sample info 
tmp_SampleInfo = readtable([Allen_dir 'ahba_dme_hcp_top8kgenes_scores.csv']);
Sample_label = tmp_SampleInfo.label;
Sample_NumLable = table2array(tmp_SampleInfo(:,1));
Gene_DME_Component = table2array(tmp_SampleInfo(:,2:4));

% read GO label and screen 
if ifloadGO
    load([Group_dir 'GO_label.mat'],'GO_categories');
else
    GO_Label = readtable([Allen_dir 'expression/GO_label.csv']);
    GO_categories = sunfun_disposalGOlabel(GO_Label);
    save([Group_dir 'GO_label.mat'],'GO_Label','GO_categories');
end

tmp_GOscreen = [];
for i = 1:size(GO_categories,1)
    if numel(GO_categories.gene_label{i}) >= 20 && strcmp(GO_categories.go_type{i},'biological_process')
        tmp_GOscreen = [tmp_GOscreen; 1];
    else
        tmp_GOscreen = [tmp_GOscreen; 0];
    end
end
GO_categories = GO_categories(tmp_GOscreen>0,:);


%% neuropiptide
neuropiptide_list = readtable([Refer_dir 'neuropeptide_gene_list.xlsx']);
tt = zeros(size(neuropiptide_list.Gene));
for i = 1:numel(GeneID)
    for j = 1:numel(tt)
        if strcmp(GeneID{i},neuropiptide_list.Gene{j})
            tt(j) = i;
        end
    end
end
tmp = table(tt,'VariableNames',{'Gene_id'});
neuropiptide_list = [neuropiptide_list tmp];
neuropiptide_list = neuropiptide_list(tt>0,:);





%% read cell type data
CT_categories1 = load([CellType_dir FileName_CellType1]).Data_CellType;
CT_categories1 = sortrows(CT_categories1,1);

CT_categories2 = load([CellType_dir FileName_CellType2]).Data_CellType;
% tmpOrder = [7 13 9 6 5 10 4 11 3 12 14 15 16 17 2 1 8];
% CT_categories2 = CT_categories2(tmpOrder,:);

CT_categories3 = load([CellType_dir FileName_CellType3]).Data_CellType;




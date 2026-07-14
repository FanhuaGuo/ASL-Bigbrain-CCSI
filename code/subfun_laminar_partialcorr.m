function [ r_pc , spin_p , explained_real , explained_null ] = subfun_laminar_partialcorr( Data1 , Data2 , Data1null , mask , CorrType )

% remove roi by mask
Data1 = Data1(mask,:);
Data2 = Data2(mask,:);
for i = 1:numel(Data1null)
    Data1null{i} = Data1null{i}(mask,:);
end

% compute real data common component
Data1 = zscore(Data1);
[Data1_common explained_real] = subfun_PCAtoPC(Data1);
Data1_common = Data1_common(:,1);

% compute partial correlation r and spin-null r
r = corr(Data1_common,Data2,'Type',CorrType);
r_pc = partialcorr(Data1,Data2,Data1_common,'Type',CorrType);
cn = {}; cn{numel(Data1null)} = [];
rn = [];
explained_null = [];
for i = 1:size(Data1null{1},2)
    tmp = [];
    for j = 1:numel(Data1null)
        tmp = [tmp Data1null{j}(:,i)];
    end
    tmp = zscore(tmp);
    [tmp_common, tExp] = subfun_PCAtoPC(tmp);
    explained_null = [explained_null; tExp(1)];
    tmp_common = tmp_common(:,1);
    tr = corr(tmp_common,Data2,'Type',CorrType);
    rn = [rn; tr];
    tt = partialcorr(tmp,Data2,tmp_common,'Type',CorrType);
    for j = 1:numel(Data1null)
        cn{j} = [cn{j}; tt(j,:)];
    end
end

% compute the p value from spin-null 
spin_p = [];
for i = 1:numel(Data1null)
    p = sum(abs(r_pc(i,:))<abs(cn{i}))/size(Data1null{1},2);
    spin_p = [spin_p; p];
end
r_pc = [r; r_pc];
spin_p = [sum(abs(r)<abs(rn))/size(Data1null{1},2); spin_p];



end


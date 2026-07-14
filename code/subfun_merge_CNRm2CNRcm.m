function [ CNRcm ] = subfun_merge_CNRm2CNRcm( data )

CMRcm_fisherz = zeros(size(data{1}));
for k = 1:numel(data)
    tmpCNRm = data{k};
    
    row_m = mean(tmpCNRm, 2, 'omitnan');
    col_m = mean(tmpCNRm, 1, 'omitnan');
    all_m = mean(tmpCNRm(:), 'omitnan');
    C0 = tmpCNRm - row_m - col_m + all_m;
    row_mean = mean(C0, 2, 'omitnan');
    row_std  = std(C0, 0, 2, 'omitnan');
    Cnorm = (C0 - row_mean) ./ (row_std + eps);
    Cnorm(1:size(Cnorm,1)+1:end) = 0;      % 对角置0

    S = corr(Cnorm', 'type','Spearman', 'rows','pairwise');
    S(~isfinite(S)) = 0;
    S(S<0) = 0;                              % 负相关置0（非负核适合 DME）

    CMRcm_fisherz = CMRcm_fisherz + subfun_fisher_z(S);
end
CMRcm_fisherz = CMRcm_fisherz/numel(data);

CNRcm = subfun_fisher_z_inverse(CMRcm_fisherz);

end
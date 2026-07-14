function [ PS1_cs ] = subfun_ProfileSimilarity_mri2mri( DataScan , DataScan2 )
%% profile similarity
[N_region, N_layer] = size(DataScan{1});

% delete the layer 1 (it's about 10% of cortical thickness)
Num_Layer1 = max([1 round(N_layer/10)]);

% compute cosine similarity
PS1_cs = [];
for k = 1:numel(DataScan)
    % Normalize each region's vector
    norm_Data = zscore(DataScan{k}(:,1:end-Num_Layer1)')';
    norm_CellDensity = zscore(DataScan2{k}(:,1:end-Num_Layer1)')';

    % 初始化结果
    cos_sim = zeros(N_region,1);   % cosine similarity

    for r = 1:N_region
        v1 = norm_Data(r, :);
        v2 = norm_CellDensity(r, :);

        % cosine similarity
        cos_sim(r) = dot(v1, v2) / (norm(v1) * norm(v2));
    end
    
    PS1_cs = [PS1_cs cos_sim];
end


end


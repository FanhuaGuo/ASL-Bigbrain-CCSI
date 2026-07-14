function [ real_matrix , perm_matrix ] = subfun_GenerateNetworkCalcMatrix_PermMatrix( Par_parcels , Net_parcels )
%
%% merge and delete 0 value
parcels = [Par_parcels Net_parcels];
parcels(parcels(:,1)==0,:) = [];
parcels(parcels(:,2)==0,:) = [];
Max_label = [max(Par_parcels),max(Net_parcels)];

%% calc real matrix
real_matrix = calculate_matrix(parcels,Max_label);

%% calc perm matrix
rng(1234);  % fix seed for repeatability
perm_matrix = {};
perm_times = 10000;
for i = 1:perm_times
    if mod(i,1000) == 0
        fprintf(['finished ' num2str(i) ' times\n']);
    end
    tmp = parcels(:,2);
    tmp = tmp(randperm(length(tmp)));
    perm_matrix{i} = calculate_matrix([parcels(:,1) tmp],Max_label);
end


end



function Out = calculate_matrix( labels , Max_label )
Out = zeros(Max_label)';
for i = 1:Max_label(2)
    tmp = labels(labels(:,2)==i,1);
    for j = 1:Max_label(1)
        Out(i,j) = sum(tmp==j);
    end
end
end

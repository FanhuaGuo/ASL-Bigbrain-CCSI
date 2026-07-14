function [ MPC ] = subfun_computeMPC( inputData1 )
% inputData1 should be ROIxlayer

cvm = corr(inputData1', 'Type', 'Pearson');
pcm = corr(inputData1',mean(inputData1,1)','Type', 'Pearson');

MPC = zeros(size(cvm));
for i = 1:size(MPC,1)
    for j = 1:size(MPC,2)
        MPC(i,j) = (cvm(i,j) - pcm(i)*pcm(j)) / sqrt((1-pcm(i)^2) * (1-pcm(j)^2));
    end
end

MPC(MPC<0) = 0;
MPC = log(MPC+1);

end


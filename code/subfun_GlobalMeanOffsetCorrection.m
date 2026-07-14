function [ outData ] = subfun_GlobalMeanOffsetCorrection( inData )
% correct left-right bias (Global Mean Offset Correction)

globalMean = nanmean(inData,1);
leftMean = nanmean(inData(1:end/2,:),1);
rightMean = nanmean(inData(end/2+1:end,:),1);

leftOffset = leftMean - globalMean;
rightOffset = rightMean - globalMean;

outData = inData;
outData(1:end/2,:) = inData(1:end/2,:) - leftOffset;
outData(end/2+1:end,:) = inData(end/2+1:end,:) - rightOffset;

end


function [ outputD ] = subfun_merge_bil2uni( inputD )

if size(inputD,2)>1
    outputD = (inputD(:,1:end/2)+inputD(:,end/2+1:end))/2;
else
    outputD = (inputD(1:end/2)+inputD(end/2+1:end))/2;
end

end


function [ OutTable ] = subfun_write_to_MapPlot_byR( MapData , DataLabel , atlasInfo , SaveName )

label = {};
score = -100*ones(size(atlasInfo,1),1);
score(DataLabel) = MapData;
for i = 1:size(atlasInfo,1)
    if atlasInfo.LR{i} == 'L'
        prefix = 'lh_L_';
    else
        prefix = 'rh_R_';
    end
    switch atlasInfo.region{i}
        case '7Pl'
            label{i} = [prefix '7PL'];
        otherwise
            label{i} = [prefix atlasInfo.region{i}];
    end
    switch score(i)
        case -100
            score(i) = 0;
        case 0
            score(i) = 0.0001;
    end
end
label{size(atlasInfo,1)+1} = 'lh_???';
label{size(atlasInfo,1)+2} = 'rh_???';
score(size(atlasInfo,1)+1) = 0;
score(size(atlasInfo,1)+2) = 0;

OutTable = table(label',score,'VariableNames',{'label','score'});
writetable(OutTable, SaveName);

end


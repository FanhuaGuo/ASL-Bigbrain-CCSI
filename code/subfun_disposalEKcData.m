function [ EKc_out , EKcw_out ] = subfun_disposalEKcData( EKc_raw )

EKc_out = [];
EKcw_out = [];
for i = 1:360
    tmp = EKc_raw(EKc_raw(:,2)==i,1);
    ttmp = [];
    for k = 1:5
        ttmp = [ttmp sum(tmp==k)];
    end
    [~,loc] = max(ttmp);
    EKc_out = [EKc_out; loc];
    EKcw_out = [EKcw_out; ttmp/sum(ttmp)];
end


end


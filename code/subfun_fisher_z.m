function [ fisher_z ] = subfun_fisher_z( r )

fisher_z = log( (1+r) ./ (1-r) )/2;

end


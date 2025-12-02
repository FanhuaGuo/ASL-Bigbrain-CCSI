function [ r ] = subfun_fisher_z_inverse( z )

r = (exp(2 .* z) - 1) ./ (exp(2 .* z) + 1);

end


function [ Output , explained ] = subfun_PCAtoPC( Input )


[coeff,score,latent,tsquared,explained,mu] = pca(Input);
Output  = Input*coeff;



end


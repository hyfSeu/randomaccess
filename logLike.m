function costOut = logLike(Y,noise_power,Zhat,Pvar)
%compute log likelihood
% For sum-product GAMP, compute
%   E( log p_{Y|Z}(y|z) ) with z = CN(phat, pvar)
% For max-sum GAMP compute
%   log p_{Y|Z}(y|z) @ z = phat
    temp = (abs(Y-Zhat).^2+Pvar)./noise_power;
    costOut = -temp;
function [ phiVal ] = phiMerton (w,param)
% phiMerton is characteristic function for Merton's 
% jump-diffusion model
    muJ=param.muJ;
    sigmaJ=param.sigmaJ;
    lambda=param.lambda;
    phiVal=exp(param.T.*(i.*param.mu.*w...
            -0.5*w.^2 .* param.sigmaBS.^2 ...
        +lambda.*(exp(i*w*muJ-0.5*sigmaJ^2.*w.^2)-1)));
end
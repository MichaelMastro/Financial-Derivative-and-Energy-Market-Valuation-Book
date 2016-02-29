function [ phiVal ] = phiBS (w,param)
% phiBS is characteristic function for Gaussian 
% log-return
    phiVal=exp(param.T.*(i.*param.mu.*w...
            -0.5*w.^2 .* param.sigmaBS.^2));
end
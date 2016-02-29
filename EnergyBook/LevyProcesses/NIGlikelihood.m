function [ like ] = NIGlikelihood (param )
% NIGlikelihood calculates negative likelihood
% of parameters of NIG process for log-return dX
% param can accept 3 (+ mu=0) or 4  values
global dX

alphaG=param(1); 
betaG=param(2); 
deltaG=param(3);
if (nargin == 4)  
    muG=param(4);
else
    muG=0;
end

n=length(dX);
z=sqrt( deltaG^2+(dX-muG).^2);
zalpha=z*alphaG;

like = n*log(alphaG*deltaG/pi)...
    + n*deltaG*sqrt(alphaG^2^2-betaG^2)...
    +sum( betaG*(dX-muG) + log( besselk(1,zalpha)./z));
like=-like;

% Crude technique to force alpha and delta
% to be greater than zero 
if (alphaG<=0), like=like+1e8; end
if (deltaG<=0), like=like+1e8; end
    
end



function [c PI1 PI2] = BlackScholesCall (K,S,T,vol,r,d);
% BlackScholesCall uses classic semi-analytical equation
    d1 = ( log(S./K)+ (r-d+vol^2/2).*T)./ (vol.*sqrt(T));
    d2 = d1-vol*sqrt(T);
    PI1 = myNormCDF(d1);
% PI1=CDF(d1) = probability of finishing in the money for
% risk neutral Martingale measure with Stock as Numeraire 
    PI2 = myNormCDF(d2);
% PI2=CDF(d2) = probability of finishing in the money for
% risk neutral Martingale measure with riskless Bond Numeraire 
    c = exp(-d*T).*S.*myNormCDF(d1)-exp(-r*T).*K.*myNormCDF(d2);
end

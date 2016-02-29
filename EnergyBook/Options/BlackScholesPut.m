function p = BlackScholesPut (K,S,T,vol,r,d);
    d1 = ( log(S./K)+ (r-d+vol^2/2).*T)./ (vol.*sqrt(T));
    d2 = d1-vol*sqrt(T);
    p=exp(-r*T).*K.*myNormCDF(-d2)-exp(-d*T).*S.*myNormCDF(-d1);
end
function ncdf = myNormCDF (x)
% ncdf=  0.5*(1+erf(x/sqrt(2)));
ncdf = 0.5*erfc(-x/sqrt(2));
end


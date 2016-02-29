function c = MertonSeriesCall(K,S,T,volBS,...
                        r, d, muJ, sigmaJ, lambda);
% MertonSeriesCall based on Series Approximation where
% n is the probability of n jumps in one time period
% lambda is the intensity parameter = mean number of jumps
% in one time period 
%k=exp(muJ+0.5*sigmaJ^2)-1 = mean relative asset jump size
    N=10; % Max Number of Possible jumps in one time period
    lambdaBar=lambda*exp(muJ+0.5*sigmaJ^2);
    c=0;
% Can write following more compactly using 
% exp(muJ+0.5*sigmaJ^2)= k+1
for n=0:N
    volM = sqrt(volBS^2+n*sigmaJ/T);
    rn = r - lambda*(exp(muJ+0.5*sigmaJ^2)-1)...
        + n*(muJ+0.5*sigmaJ^2)/T;
    c = c + ((lambdaBar*T)^n / factorial(n))...
        * BlackScholesCall (K,S,T,volM,rn,0);
end
    c = c*exp(-lambdaBar*T);
end
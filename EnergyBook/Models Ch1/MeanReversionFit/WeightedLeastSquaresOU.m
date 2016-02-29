function [mu, lambda, sigma,  CalcSlope, CalcIntercept, StdDev]...
    = WeightedLeastSquaresOU (S,delta,sigma)
% WeightedLeastSquaresOU performs a weighted or unweighted 
% least squares fit to Orstein Uhlenbeck process
% Derivation in Press, Flannery, Teukolsky, Vetterling, 
% Numerical Recipes
% as well as an unweighted estimation procedure 
% by M.A. van den Berg available at www.sitmo.com

x=S(1:end-1);
y=S(2:end);

n= length (y);

if nargin < 3, sigma = ones(1,n); end

S   = sum (1./sigma.^2);
Sx  = sum (x./sigma.^2);
Sy  = sum (y./sigma.^2);
Sxx = sum (x.^2./sigma.^2);
Sxy = sum (x.*y./sigma.^2);
Syy = sum (y.^2./sigma.^2);

CalcSlope    =(S.*Sxy-Sx.*Sy)  /(S.*Sxx-Sx.^2);
CalcIntercept=(Sxx.*Sy-Sx.*Sxy)/(S.*Sxx-Sx.^2);
StdDev=sqrt((S*Syy - Sy.^2 - (CalcSlope.*(S.*Sxy-Sx.*Sy)))...
    / (S*(S-2)) ); 
%%UnWeighted
%=sqrt((n*Syy - Sy.^2 - (CalcSlope.*(n.*Sxy-Sx.*Sy))) / (n*(n-2)) ); 

lambda = -log(CalcSlope)/delta;
mu = CalcIntercept/(1-CalcSlope);
sigma = StdDev * sqrt(-2*log(CalcSlope)/(delta*(1-CalcSlope^2)));

end

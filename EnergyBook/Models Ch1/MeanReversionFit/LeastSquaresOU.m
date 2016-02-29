function [mu, sigma, lambda] = LeastSquaresOU (S,delta)

%[a, b, StandardDeviation] = LeastSquaresLine (S(1:end-1),S(2:end))

x=S(1:end-1);
y=S(2:end);

n= length (y);
xMean=mean(x);
yMean=mean(y);

SxSy= sum (x.*y) - n*xMean*yMean;
SySy= sum (y.*y) - n*yMean*yMean;
SxSx= sum (x.*x) - n*xMean*xMean;

CalcSlope=SxSy/SxSx;
CalcIntercept=yMean-CalcSlope*xMean;
StdDev=sqrt((SySy-CalcSlope*SxSy)/(n-2)); %=sqrt((SySy-(SxSy^2/SxSx))/(n-2))

fprintf(1, ' LS Calc.\t %6.2f Slope \t %6.2f Intercept %6.2f Standard Deviation  \n',...
    CalcSlope, CalcIntercept, StdDev);

lambda = -log(CalcSlope)/delta;
mu = CalcIntercept/(1-CalcSlope);
sigma = StdDev * sqrt(-2*log(CalcSlope)/(delta*(1-CalcSlope^2)));

end

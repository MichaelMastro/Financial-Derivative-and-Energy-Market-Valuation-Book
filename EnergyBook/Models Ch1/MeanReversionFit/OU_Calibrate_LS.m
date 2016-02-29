function [mu,sigma,lambda] = OU_Calibrate_LS(S,delta)
  n = length(S)-1;

  Sx  = sum( S(1:end-1) );
  Sy  = sum( S(2:end) );
  Sxx = sum( S(1:end-1).^2 );
  Sxy = sum( S(1:end-1).*S(2:end) );
  Syy = sum( S(2:end).^2 );

  a  = ( n*Sxy - Sx*Sy ) / ( n*Sxx -Sx^2 );
  b  = ( Sy - a*Sx ) / n;
  sd = sqrt( (n*Syy - Sy^2 - a*(n*Sxy - Sx*Sy) )/n/(n-2) );

  lambda = -log(a)/delta;
  mu     = b/(1-a);
  sigma  =  sd * sqrt( -2*log(a)/delta/(1-a^2) );
end

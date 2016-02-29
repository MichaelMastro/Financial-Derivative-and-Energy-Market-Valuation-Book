function [mu,lambda, sigma] = weightedML(S,delta,sigma)
%weightedML performs a maximum likelihood estimation of the 
%the parameters in a Orstein Uhlenbeck process 
%Also see detailed descriptions of the 
%unweighted estimation  
%by M.A. van den Berg available at www.sitmo.com
%and Weisstein, Eric W. MathWorld, Wolfram Research,
%http://mathworld.wolfram.com

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

%In derivation, two equations and two unknowns are available 
%for mu and lambda. 
%Sigma is directly solvable once mu and lambda are calculated
  mu  = (Sy*Sxx - Sx*Sxy) / ( S*(Sxx - Sxy) - (Sx^2 - Sx*Sy) );
  lambda = -log( (Sxy - mu*Sx - mu*Sy + S*mu^2) /...
      (Sxx -2*mu*Sx + S*mu^2) ) / delta;
  a = exp(-lambda*delta);
  sigmah2 = (Syy - 2*a*Sxy + a^2*Sxx - 2*mu*(1-a)*(Sy - a*Sx)...
      + S*mu^2*(1-a)^2)/S;
  sigma = sqrt(sigmah2*2*lambda/(1-a^2));
end

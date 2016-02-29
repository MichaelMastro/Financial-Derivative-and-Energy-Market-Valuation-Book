function [ dIG ] = InvGaussianGen ( a,b )
% GammaRand generates Gaussian Random variables 
% IG( at, b)
% Pseudo-Code and description in Schoutens
% Levy Processes in Finance
% of Michael, Schucany, and Haas generator

if (nargin <2), b=49; end
if (nargin <1), a=1; end

v=randn;
y=v.^2;
x=(a/b)+y/(2*b^2)-sqrt(4*a*b*y+y^2)/(2*b^2);
u=randn;
if (u <= (a/(a+x*b)))
    dIG=x;
else
    dIG=a^2/(b^2*x);
end

end


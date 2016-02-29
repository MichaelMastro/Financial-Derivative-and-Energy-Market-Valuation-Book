function [ dX, dIG ] = NormInvGaussianGen(alpha, beta, delta, dt )
% NormInvGaussianGen calculates dX of a Inverse-Gaussian 
% time process subordinated Brownian motion
if (nargin < 4), dt = 1; end
if (nargin < 3), theta = 10; end
if (nargin < 2), alpha = 50; end
if (nargin < 1), beta = -10; end
   
a=1; b=delta*sqrt(alpha^2-beta^2);
dIG=InvGaussianGen (a.*dt,b); % subordinate Inv Gaussian time
%subordinated drift - diffusion
dX=beta.*delta.^2.*dIG + delta*sqrt(dIG).*randn; 
end

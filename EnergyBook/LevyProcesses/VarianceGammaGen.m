function [ dX, dG ] = VarianceGammaGen( sigma, v, theta, dt )
% VarianceGammaGen calculates dX of a Variance Gamma 
% subordinated Brownian motion
if (nargin < 4), dt =1; end
if (nargin < 3), theta =1; end
if (nargin < 2), v =1; end
if (nargin < 1), sigma =1; end
   
a=1/v; b=1/v;
dG=GammaRand (a.*dt,b); % subordinate Gamma time process
dX=theta.*dG+sigma*sqrt(dG).*randn; %subordinated drift - diffusion
end

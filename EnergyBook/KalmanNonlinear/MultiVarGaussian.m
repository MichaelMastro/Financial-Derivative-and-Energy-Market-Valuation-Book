function [ ] = MultiVarGaussian( ) 
% MultiVarGaussian Plots 2D Gaussian Distribution and 
% Gauss-Hermite Quadrature Points. A stochastic decoupling
% technique based on Cholesky Decomposition and SVD is used to 
% eliminate correlation
% Matrix points are rotated via Cholesky Decomposition:
% [R D] = chol(Sigma); L=R'; xnew=L*x;
% or SVD: [U,S,V] = svd(Sigma); xnew=U*sqrt(S)*x; 
close all
clc

set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold') 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

rho = 0.6
sig= [1 rho ;rho 1]
mu = [0;0]

MaxRange=3;
MaxNum=40
M=linspace(-MaxRange,MaxRange,MaxNum);
N=linspace(-MaxRange,MaxRange,MaxNum);
[X,Y] = meshgrid(M,N);

for m = 1:MaxNum
    for n=1:MaxNum
        x=[M(m); N(n)];
        densX(m,n)  = MultiPDF( x, mu, sig ); 
    end
end

% Calculate Gauss-Hermite Weights and Quadrature Points
NumP=5
[xp, wp ] = HermiteWeightAndRoots (NumP, 0)
% Sort Quadrature Points
[xp ind] = sort(xp); % [Xp,Yp] = meshgrid(xp,xp);
% Sort Weights based on Quadrature Points
tempwp=wp;
for k =1:length(wp)
    wp(k)=tempwp(ind(k));
end
W=wp*wp'; % 2D grid of weights

% Cholesky Decomposition
[R D] = chol (sig)
L=R'

figure
subplot (1,2,1)
contour(X,Y,densX,40)
hold on
for p1=1:NumP
    for p2=1:NumP
        x=[xp(p1) xp(p2)];
        xRot = L*x';
        plot(xRot(1),xRot(2),'o','MarkerSize',2+20*(W(p1,p2))) 
        hold on
    end
end
axis ([-MaxRange-1, MaxRange+1,-MaxRange-1,MaxRange+1])
title('Cholesky Stochastic Decoupling')
txtstr=['\rho =' num2str(rho)];
text(2,-2,txtstr);
hold off

[U,S,V] = svd(sig)
A=U*sqrt(S)
subplot (1,2,2)
contour(X,Y,densX,40)
hold on
for p1=1:NumP
    for p2=1:NumP
        x=[xp(p1) xp(p2)];
        xRot = A*x';
        plot(xRot(1),xRot(2),'o','MarkerSize',2+20*(W(p1,p2))) 
        hold on
    end
end
axis ([-MaxRange-1, MaxRange+1,-MaxRange-1,MaxRange+1])
title('SVD Stochastic Decoupling')
text(2,-2,txtstr)
hold off
end

function [ Dens ] = MultiPDF( x, mu, sig ) 
% Multivariate Probability Density Function
    k=rank (sig);
    Dens = ( exp(-0.5*(x-mu)'*inv(sig)*(x-mu) ) )/...
        ( (2*pi)^(k/2)*sqrt(det(sig)));
end

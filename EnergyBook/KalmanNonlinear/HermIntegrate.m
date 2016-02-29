function [ ] = HermIntegrate(  )
% HermIntegrate implements Gauss-Hermite Quadrature for
% Definite (-inf to +inf) Integral ( f(x)*exp(-x^2) ) 
%       = sum ( w_i*f(x_i)*h(x_i) )
% for three separate functions f(x) = x, x^2, e^2x
% x_i = Gauss-Hermite quadrature points 
% w_i = Gauss-Hermite quadrature weights 
% Analytical Definite Integrals (-inf to +inf) are available for
% xexp(-x^2), x^2exp(-x^2), e^2x*exp(-x^2)
% Quadrature Value compared to Analytical Value as a function of
% number of quadrature points, NumP.

set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold') 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);
close all
clear all
clc

if (nargin == 0), MaxNumP=10; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp ('Indefinite Integral (x*exp{-x^2}) = -(1/2)*exp{-x^2}')
disp (' Definite Integral from -inf to inf: (x*exp{-x^2}) = 0')  
Analytical= 0
% MyInf=1e1/eps; % use as approx inf for Matlab quadl or quad
% F = @(x) x.*exp(-x.^2);
% Ql = quadl(F,-MyInf,MyInf) 
Diff1=NaN*zeros(1,MaxNumP)
for NumP=2:MaxNumP;
    [xp, wp ] = HermiteWeightAndRoots (NumP, 0);
    Int=sum(wp.*(xp));
    Diff1(NumP)=Analytical-Int;
    fprintf (1,...
    'n = %g, G-H Quad x*exp{-x^2} = %e; Error = %e \n',...
        NumP, Int, Diff1(NumP))
end
fprintf (1,'\n\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Diff2=NaN*zeros(1,MaxNumP);
%disp ('Indefinite Integral (x^2*exp{-x^2}) = exp(x)*(x^2-2x+2)')
fprintf (1, ' Definite Integral from -inf to inf:')
fprintf (1, '  x^2*exp{-x^2} = 0.5*sqrt(pi)\n')
Analytical= 0.5*sqrt(pi)
% F = @(x) x.^2.*exp(-x.^2);
%Ql = quadl(F,-MyInf,MyInf) 
for NumP=2:MaxNumP;
    [xp, wp ] = HermiteWeightAndRoots (NumP, 0);
    Int=sum(wp.*(xp.^2));
    Diff2(NumP)=Analytical-Int;
    fprintf (1,...
    'n = %g, G-H Quad x^2*exp{-x^2} = %e; Error = %e \n',...
        NumP, Int, Diff2(NumP))
end
fprintf (1,'\n\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf (1,' Definite Integral from -inf to inf:')
fprintf (1,' exp(-2bx)*exp{-x^2} = exp(b^2)*sqrt(pi)\n')
Analytical= exp(1^2)*sqrt(pi) % b=1
% F = @(x) exp(-2*x).*exp(-x.^2);
% Ql = quadl(F,-MyInf,MyInf) 
Diff3=NaN*zeros(1,MaxNumP);
for NumP=2:MaxNumP;
    [xp, wp ] = HermiteWeightAndRoots (NumP, 0);
    Int=sum(wp.*(exp(-2*xp)));
    Diff3(NumP)=Analytical-Int;
    fprintf (1,...
    'n = %g, G-H Quad exp(-2bx)*exp{-x^2} = %e; Error = %e \n',...
        NumP, Int, Diff3(NumP))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Error in Quadrature Relative to Analytical Value
figure
subplot(3,2,1); 
plot (1:MaxNumP,Diff1);
xlabel ('n = Quadrature Points'); ylabel ('Error')
title('\int_{-\infty}^{\infty}xe^{-x^2} = \Sigmax_ih(x_i)')
subplot(3,2,2); 
semilogy (1:MaxNumP,abs(Diff1));
xlabel ('n = Quadrature Points'); ylabel ('Abs. Error')
title('\int_{-\infty}^{\infty}xe^{-x^2} = \Sigmax_ih(x_i)')

subplot(3,2,3);
plot (1:MaxNumP,Diff2);
xlabel ('n = Quadrature Points'); ylabel ('Error')
title('\int_{-\infty}^{\infty}x^2e^{-x^2} = \Sigmax_i^2h(x_i)')
subplot(3,2,4); 
semilogy (1:MaxNumP,abs(Diff2));
xlabel ('n = Quadrature Points'); ylabel ('Abs. Error')
title('\int_{-\infty}^{\infty}x^2e^{-x^2} = \Sigmax_i^2h(x_i)')

subplot(3,2,5); 
plot (1:MaxNumP,Diff3);
xlabel ('n = Quadrature Points'); ylabel ('Error')
title(...
'\int_{-\infty}^{\infty}e^{-2x}e^{-x^2} = \Sigmae^{-2x_i} h(x_i)')
subplot(3,2,6); 
semilogy (1:MaxNumP,abs(Diff3));
xlabel ('n = Quadrature Points'); ylabel ('Abs. Error')
title(...
 '\int_{-\infty}^{\infty}e^{-2x}e^{-x^2} = \Sigmae^{-2x_i} h(x_i)')

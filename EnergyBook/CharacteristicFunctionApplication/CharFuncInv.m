function pdf = CharFuncInv (); 
%CharFuncInv produces PDF via Fourier Inversion of
%alpha-stable characteristic function
close all; clear all; clc
set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold') 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

p.t = 1/252;
p.mu=0.1; %or risk neutral pricing use r-0.5*sigma^2;
p.alpha=1.9;
p.beta=0.1;
p.sigma =0.3;

%characteristic function discretization
%Increase Max Omega for fatter tails
%Certainly there is a better approach for OmegaEnd
OmegaEnd=round(2/(p.sigma*p.t*(p.alpha^2)));
%frequency grid size for char func
omega=linspace(0,OmegaEnd,140)';
DeltaOmega=omega(2)-omega(1)    

x = linspace (-0.2,0.2,91);
Nx=length (x)
pdf=zeros(Nx,1);

phin=phi_alpha(omega,p);
%Show PDF of one point (x_n) is integral of 
%Real component of characteristic function times
%e^(-i x_n omega) 
%where e^(-i x_n omega) = cos() + i sin () 
figure
subplot (3,2,1); plot(omega, real(phin));
title('Real Char. Funct.'); xlabel('\omega'); 
ylabel('Real(\phi_\omega)'); axis tight;

subplot (3,2,2); plot(omega, imag(phin));
title('Imag Char. Func.'); xlabel('\omega'); 
ylabel('Imag(\phi_\omega)'); axis tight;

subplot (3,2,3); plot(omega, real(exp(-i*x(1)*omega)));
title('Real x_1 Oscillation'); xlabel('\omega');
ylabel('Real(e^{-ix_1\omega})'); axis tight;

subplot (3,2,4); plot(omega, imag(exp(-i*x(1)*omega)));
title('Imaginary x_1 Oscillation'); xlabel('\omega');
ylabel('Imag(e^{-ix_1\omega})'); axis tight;

subplot (3,2,5:6); 
area(omega, real(exp(-i*x(1)*omega).*phin));
title('PDF(x_1) ~ Area '); xlabel('\omega');
ylabel('Real(\phi(\omega)e^{-ix_1\omega})');

%show similar area (integral) graphs for several x_n
figure
for k=1:Nx
    y=exp(-i*x(k)*omega).*phin;
    y=real(y);
    if (k==1)
        title('PDF(x_n) ~ Area ');
        subplot(4,1,1); area(omega,y);
        numb2=num2str(x(k)); numb=num2str(k);
        ylabel(['x_{' numb '}=' numb2]);     
    end
    if (k == round(Nx/6) )
        subplot(4,1,2); area(omega,y);
        numb2=num2str(x(k),2); numb=num2str(k);
        ylabel(['x_{' numb '}=' numb2]);        
    end
    if (k == round(Nx/3) )
        subplot(4,1,3); area(omega,y);
        numb2=num2str(x(k),2); numb=num2str(k);
        ylabel(['x_{' numb '}=' numb2]);        
    end
    if (k == round(Nx/2) )
        subplot(4,1,4); area(omega,y);
        numb2=num2str(x(k)); numb=num2str(k);
        ylabel(['x_{' numb '}=' numb2]);      
    end
    %integrate area under (or below) curve
    I=sum(y)-0.5*(y(1)+y(end)); % trapezoidal rule
    pdf(k)=I*DeltaOmega/pi; %output
end
xlabel('\omega'); 

%plot pdf, can display numerical oscillations if
%omega_end too small, that is, cut off 
%high frequency oscillations 
figure
subplot(2,1,1); plot (x,pdf)
xlabel('x'); ylabel('PDF(x)');
subplot(2,1,2); semilogy (x,abs(pdf))
xlabel('x'); ylabel('PDF(x)');
end

%



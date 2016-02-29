% cdfCharFuncFFTscript produces CDF via Fourier Inversion 
% of alpha-stable characteristic function
close all; clear all; clc
set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold') 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);


p.t = 1/12;
p.mu=0;
p.alpha=2;
p.beta=0;
p.sigma =0.5;


N=2^7;
OmegaEnd=50
figure
maxLoop=3
for loop=1:maxLoop
    eta=0.1*(loop)^3
    [CDF FFTinput omega x] = cdfCharFuncFFT(p,N,OmegaEnd,eta);
    
    ind1=loop;
    subplot (3,maxLoop,ind1); 
    ETAtitle=num2str(eta); 
    ETAtitle=['\eta=' ETAtitle ]; 
    
    plot(omega, real(FFTinput)); title(ETAtitle); 
    ylabel('Re(FFTin)'); xlabel('\omega');  axis tight
    
    ind2=loop+maxLoop;
    subplot (3,maxLoop,ind2); 
    plot(omega, imag(FFTinput))
    ylabel('Im(FFTin)'); xlabel('\omega');  axis tight
    
    ind3=loop+maxLoop*2;
    subplot (3,maxLoop,ind3);
    plot (x, real(CDF)) %semilogy
    ylabel('CDF(x)'); xlabel('x'); % axis tight
    ylim([0.01 0.99])
end




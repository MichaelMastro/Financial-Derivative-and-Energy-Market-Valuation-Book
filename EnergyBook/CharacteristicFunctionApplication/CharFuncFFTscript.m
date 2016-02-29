% CharFuncFFTscript  
%CharFuncInv produces PDF via Fourier Inversion of
%alpha-stable characteristic function
close all; clear all; clc
set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold') 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

p.t = 1/12;
p.mu=0.1;
p.alpha=1.7;
p.beta=0.4;
p.sigma =1;


%%%%%%% Compare change in number of data points
OmegaEnd=30
figure
maxLoop=3
for loop=1:maxLoop
    N=2^(6+loop)
    [PDF FFTinput omega x] = CharFuncFFT(p,N,OmegaEnd);
    
    ind1=loop;
    subplot (3,maxLoop,ind1); 
    Ntitle=num2str(N); OmegaEndTitle=num2str(OmegaEnd);
    Ntitle=['Max \omega=' OmegaEndTitle '; N=' Ntitle];
    plot(omega, real(FFTinput)); title(Ntitle); 
    ylabel('Re(FFTin)'); xlabel('\omega');  axis tight
    
    ind2=loop+maxLoop;
    subplot (3,maxLoop,ind2); 
    plot(omega, imag(FFTinput))
    ylabel('Im(FFTin)'); xlabel('\omega');  axis tight
    
    ind3=loop+maxLoop*2;
    subplot (3,maxLoop,ind3);
    semilogy (x, real(PDF))
    ylabel('PDF(x)'); xlabel('x'); axis tight
end

%%%%%% compare change in omega max  
figure
N=N/2
maxLoop=3
OmegaMid=OmegaEnd/2
for loop=1:maxLoop
    OmegaEnd=OmegaMid*loop
    [PDF FFTinput omega x] = CharFuncFFT(p,N,OmegaEnd);
    
    ind1=loop;
    subplot (3,maxLoop,ind1); 
    Wtitle=num2str(OmegaEnd); Ntitle=num2str(N)
    Wtitle=['Max \omega=' Wtitle '; N=' Ntitle];
    plot(omega, real(FFTinput)); title(Wtitle); 
    ylabel('Re(FFTin)'); xlabel('\omega');  axis tight
    
    ind2=loop+maxLoop;
    subplot (3,maxLoop,ind2); 
    plot(omega, imag(FFTinput))
    ylabel('Im(FFTin)'); xlabel('\omega');  axis tight
    
    ind3=loop+maxLoop*2;
    subplot (3,maxLoop,ind3);
    semilogy (x, real(PDF))
    ylabel('PDF(x)'); xlabel('x');  axis tight
end

%%%%%%%%% Vary alpha in fractional FFT
figure
OmegaEnd=30
N=512
maxLoop=3
for loop=1:maxLoop
    alpha=1/(loop*N)
    [PDF FFTinput omega x] = CharFuncFrFFT(p,N,OmegaEnd, alpha);
    
    ind1=loop;
    subplot (3,maxLoop,ind1); 
    Looptitle=num2str(loop); Ntitle=num2str(N)
    Looptitle=['A=1/' Looptitle 'N ; N=' Ntitle];
    %Looptitle=['\alpha=1/' Looptitle 'N ; N=' Ntitle];
    plot(omega, real(FFTinput)); title(Looptitle); 
    ylabel('Re(FFTin)'); xlabel('\omega');  axis tight
    
    ind2=loop+maxLoop;
    subplot (3,maxLoop,ind2); 
    plot(omega, imag(FFTinput))
    ylabel('Im(FFTin)'); xlabel('\omega');  axis tight
    
    ind3=loop+maxLoop*2;
    subplot (3,maxLoop,ind3);
    semilogy (x, real(PDF))
    ylabel('PDF(x)'); xlabel('x');  axis tight
end
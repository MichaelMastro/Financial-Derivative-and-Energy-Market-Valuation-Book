function [] = padFFT()
% padFFT shows the benefit of zero padding 
% the FFT input array to resolve the frequency amplitudes
% padding does not improve the frequency resolution
% but padding does decrease the frequency spacing in the
% FFT output. The example shows that the padding allows     
% the two separate peaks to be resolved

close all
clear all
clc

set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold') 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

xlength=2^6
x=1:xlength;
f1=2^3
f2=1.1*f1
y=sin(f1*2*pi*x/xlength)+sin(f2*2*pi*x/xlength);
subplot (2,2,1); plot (x,y)
xlabel('Real-Space x'); ylabel('Real-Space y'); 
axis tight
H=fft(y);
PowerH=H.*conj(H)/xlength;
subplot (2,2,2); plot((0:xlength/2),PowerH(1:xlength/2+1));
axis tight
ylabel('Power'); xlabel('Frequency');

mult=2
yp=[y zeros(1,(mult-1)*xlength)];
xp=1:(xlength*mult);
subplot (2,2,3); plot (xp/mult,yp); 
xlabel('Real-Space x'); ylabel('Zero-Padded y'); axis tight
Hp=fft(yp);
PowerHp=Hp.*conj(Hp)/xlength;
subplot (2,2,4); 
plot((0:xlength*mult/2)/2,PowerHp(1:xlength*mult/2+1));
axis tight
ylabel('Zero-Padded Power'); xlabel('Frequency');

end

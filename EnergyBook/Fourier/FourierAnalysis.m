function [ ] = FourierAnalysis( flag )
%FourierAnalysis shows DFT and FFT are equivalent but
%FFT in much faster. FFT: O(nlog2n); DFT O(n^2)
%An inverse transform is performed to show original data
%can be recovered 
%DFT: H_n= sum (h_n * exp ( -2pi*i*n*k/N) )

close all
clear all
clc

set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold') 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

%defaults to sine function for data
if (nargin == 1), %Check for Stock Data Input Flag
    dat=load('XOMprice.dat')'; %read directly 
    S=log(dat(1:end-1))-log(dat(2:end));
    SampleRate=1; %1 sample per day
    T = length(S)-1; %Format should have 1 column: price
    n=0:1/SampleRate:T;
else
    T=50;  %number of units, e.g., seconds
    SampleRate=2; % e.g., samples taken per unit (seconds)
    n = ((1:T*SampleRate)-1)/SampleRate;

    MainFreq=0.1; %=fractions of oscillations per unit(seconds)
    %angular frequency = 2*pi*MainFreq
    S=sin( 2.*pi.*n*(MainFreq));
end 

h=S; %Just so can use put into h<->H format
N = length(h);
H=zeros(1,N);
counter0=[1:N]-1;

disp('Discrete Fourier Transform');
tic
for k = 1:N
    %H(k)=0;
   % for n = 1:N %sum over all data points
    %    H(k)=H(k)+ h(n).*exp(-i.*2.*pi.*(n-1).*(k-1)./N);
    %end 
    % Partially vectorize to improve speed
    H(k)=sum(h(1:N).*exp(-1.*i.*2.*pi.*counter0.*(k-1)./N));       
end
toc

%Compare Real and Imaginary Input; Real and Imaginary DFT
%and DFT Power and Phase;
figure
freq=(0:N-1)/(N/SampleRate);
freqshift=(-N/2:1:N/2-1)/(N/SampleRate);

subplot (3,2,1); plot (n,S); axis tight;
if (nargin == 1)
    title('XOM Daily Return Input')
    ylabel('ln(S_n)-ln(S_{n-1})');
    xlabel('Time [Days]');
else
    title('Sinusoidal Real Input');
    xlabel('Unit (Time)');
end

subplot (3,2,2); plot (freq,real(H))
title('Real Fourier');
ylabel('Real(DFT)');  axis tight; 
if (nargin == 1)
    xlabel('Freq. [Cyles/Day]');
else
    xlabel('Freq. [Cyles/Unit(Time)]');
end


subplot (3,2,3); plot (freqshift,fftshift(real(H))) 
title('Shift Real(DFT)');
ylabel('fftshift');  axis tight; xlabel('Frequency');

subplot (3,2,4); plot (freq, imag(H)) 
title('Imaginary Fourier');
ylabel('Imaginary(DFT)');  axis tight; xlabel('Frequency');


subplot (3,2,5); plot (freq, abs(H)) % or H.*conj(H)/N 
title('Magnitude(DFT)');  axis tight; xlabel('frequency');
ylabel('|DFT|'); %also plot power = H*conj(H)/N
subplot (3,2,6); plot (freq, unwrap(angle(H)))   
title('Phase(DFT)'); axis tight; xlabel('frequency');
ylabel('Radians');

%{
%Show output of DFT and FFT are equivalent
%Calculate that FFT is much fast
figure
%polar form of complex number: z=rexp(i*theta)
subplot (2,2,1); plot (freq,abs(H) ) %or H.*conj(H)/N 
title('Discrete Fourier Transform');
ylabel('|FFT|');  axis tight; xlabel('frequency');
subplot (2,2,2); plot (freq, unwrap(angle(H))) 
ylabel('Phase'); axis tight; xlabel('frequency');
h=zeros(1,N);%clear h; %initialize

disp('Fast Fourier Transform');
tic
fastF=fft(S);
toc

subplot (2,2,3); plot (freq, abs(fastF)) %~abs(H)
title('Fast Fourier Transform');
ylabel('|FFT|');  axis tight; xlabel('frequency');
subplot (2,2,4); plot (freq, unwrap(angle(fastF))) 
ylabel('Phase'); axis tight; xlabel('frequency');

disp('Inverse Discrete Fourier Transform');
tic 
for n = 1:N
    %h(n)=0;
    %for k = 1:N %sum over all data points
     %   h(n)=h(n)+ H(k).*exp(i.*2.*pi.*(n-1).*(k-1)./N);
    %end
    %h(n)=h(n)/N;  % Product of iDFT and DFT must equal 1/N
    % Partially vectorize to improve speed
    invf(n)=sum ( H(1:N).*exp(i.*2.*pi.*(n-1).*(counter0)./N) ) / N;
end
toc

disp('Inverse Fast Fourier Transform');
tic
invfastF=ifft(fastF);
toc

%Show FT and inverse FT preserves data
figure
subplot (3,1,1); plot (S)
title('Original XOM'); axis tight;
if (nargin == 1)
    title('XOM')
end
axis tight; 
subplot (3,1,2); plot (real(invf))
title('DFT H^-^1(H(h))'); axis tight; 
if (nargin==1)
    ylabel('daily return');
    xlabel('time [days]');
end
subplot (3,1,3); plot (invfastF)
title('FFT H^-^1(H(h))'); axis tight;
%}
end



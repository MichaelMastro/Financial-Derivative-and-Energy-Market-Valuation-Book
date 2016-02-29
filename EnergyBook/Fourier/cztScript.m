% cztScript Calls mczt = my Chirp z-Transform function
% The mczt and Matlab's czt Routine use Bluestein and
% the Rabiner/Schafer routine based on one iFFT and 
% two FFTs. For comparison mczt can run the slower direct
% Fourier transform or a convolution. 

% The first part of the script confirms that the output is 
% identical and compares the speed of each technique.

% The second half of the script shows that the number of 
% output samples M can be more, equal to, or less than the
% the number of input data points N. Also the range of 
% frequencies can be narrowed to a region of interest via 
% the W and A parameters.

% The third (last) part of the script calls the fractional
% FFT function frfft, which is another form of the
% chirp z transform. The impact of a change in alpha (fraction)
% is shown in the angular frequency end and the resultant
% power-frequency spectrum
close all
clear all
clc

set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold') 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

N=2^9
n=(0:N-1);
% Generate Data X(n) with 3 sinusoidals
f1=N/4
f2=1.01*f1
f3=1.03*f1
X=sin(f1*2*pi*n/N)+sin(f2*2*pi*n/N)+sin(f3*2*pi*n/N);

figure
plot(n,X)
xlabel('Real-Space n'); Ylabel('Real-Space X');
title ('3 Closely Spaced Sinusoidals');

% M = number of output data points
M = N/2; %M can be <,=,or > N

fstart=N/5
fstop=N/3

frange=N;
phi=(fstop-fstart)/frange/M;
theta=fstart/frange;
W=exp(-i*2*pi*phi); % Angular Spacing
 
A=exp(i*2*pi*theta);  % Angular Starting Point

%%%%%%% Compare Computational Methods
figure %%%%%%% 

display ('Matlab FFT')
tic; y = fft (X); toc;
y=y.*conj(y)/N; 
outputM = length(y) 
k=0: (outputM-1);
subplot (4,1,1); plot (k,y); title ('FFT')
xlabel ('Frequency'); ylabel ('Power');

method ='direct'; display ('Direct FT')
tic; y = mczt ( X, M , W, A, method); toc;
outputM = length(y)
y=y.*conj(y)/N;
k=linspace(fstart,fstop, outputM); 
subplot (4,1,2); plot (k,y); title ('Direct')
xlabel ('Frequency'); ylabel ('Power');

method ='convolution';  display ('Convolution')
tic; y = mczt ( X, M , W, A, method); toc;
outputM = length(y) 
k=linspace(fstart,fstop, outputM);
y=y.*conj(y)/N;
subplot (4,1,3); plot (k,y); title ('Conv.')
xlabel ('Frequency'); ylabel ('Power');

method ='Bluestein'; display ('Bluestein')
tic; y = mczt ( X, M , W, A, method); toc;
outputM = length(y) 
k=linspace(fstart,fstop, outputM);
y=y.*conj(y)/N;
subplot (4,1,4); plot (k,y); title ('Bluestein')
xlabel ('Frequency'); ylabel ('Power');

%%% show CZT sample variation range%%%%%%%%%%%%%%%%
display ('Vary Number of Output and Input Bins')
N=2^5;
n=(0:N-1);
f1=N/4
f2=1.1*f1
X=sin(f1*2*pi*n/N)+sin(f2*2*pi*n/N);
figure
plot(n,X)
xlabel('Real-Space n'); Ylabel('Real-Space X');
title ('2 Closely Spaced Sinusoidals');

figure 
%%%%%%%%% CZT replicates FFT 
M=N;
fstart=0;
fstop=N;
frange=N;
phi=(fstop-fstart)/frange/M;
theta=fstart/frange;
W=exp(-i*2*pi*phi);
A=exp(i*2*pi*theta);

k=0:(M-1);
z=A.*W.^(-k);
subplot (3,2,1); plot (z,'*'); 
xlim([-1 1]); ylim([-1 1]);
xlabel('Real(z_n)'); ylabel('Imaginary(z_n)')
title ('N=2^5; M=2^5')

method ='Bluestein'; 
display ('Bluestein Replicates FFT')
tic; y = mczt ( X, M , W, A, method); toc;
outputM = length(y) 
k=linspace(fstart,fstop, outputM);
y=y.*conj(y)/N;
subplot (3,2,2); plot (k,y); 
title ('Bluestein Repl. FFT'); xlabel ('Frequency');


%%%%%%%%% Increase Density of points by increasing M %%%
M=4*N;
fstart=0;
fstop=N;
frange=N;
phi=(fstop-fstart)/frange/M;
theta=fstart/frange;
W=exp(-i*2*pi*phi);
A=exp(i*2*pi*theta);

k=0:(M-1);
z=A.*W.^(-k);
subplot (3,2,3); plot (z,'*'); 
xlim([-1 1]); ylim([-1 1]);
xlabel('Real(z_n)'); ylabel('Imaginary(z_n)')
title ('N=2^5; M=2^7')

method ='Bluestein'; 
display ('Increase Number of Output Bins')
tic; y = mczt ( X, M , W, A, method); toc;
outputM = length(y) 
k=linspace(fstart,fstop, outputM);
y=y.*conj(y)/N;
subplot (3,2,4); plot (k,y); 
title ('Inc. Output'); xlabel ('Frequency');

%%%%%%%%% Decrease Frequency Range
M=N;
fstart=N/3;
fstop=N/5;
frange=N;
phi=(fstop-fstart)/frange/M;
theta=fstart/frange;
W=exp(-i*2*pi*phi);
A=exp(i*2*pi*theta);

k=0:(M-1);
z=A.*W.^(-k);
subplot (3,2,5); plot (z,'*'); 
xlim([-1 1]); ylim([-1 1]);
xlabel('Real(z_n)'); ylabel('Imaginary(z_n)')
title ('N=2^5; M=2^5')

method ='Bluestein'; 
display ('Decrease Frquency Range')
tic; y = mczt ( X, M , W, A, method); toc;
outputM = length(y) 
k=linspace(fstart,fstop, outputM);
y=y.*conj(y)/N;
subplot (3,2,6); plot (k,y); 
title ('Dec. Freq Range'); xlabel ('Frequency');

figure
%%%%%%%%% fractional FFT %%%%%%%%%%%%%%%%%%%%
% Calls frfft which then calls mczt
% factional FFT has same number of input (N) and output
% bins (M); alpha parameter changes spacing in W
% and thus changes end frequency
% Starting point is always at zero frequency, i.e., A=1
alpha = 1/N;
M=N;
W=exp(-i*2*pi*alpha);
A=1;

k=0:(N-1);
z=A.*W.^(-k);
subplot (3,2,1); plot (z,'*'); 
xlim([-1 1]); ylim([-1 1]);
xlabel('Real(z_n)'); ylabel('Imaginary(z_n)')
title ('\alpha = 1/N, N=M=2^5')

display ('Fractional FFT: alpha = 1/N')
tic; y = frfft ( X, alpha); toc;
outputM = length(y) 
k=linspace(0,N^2*alpha, outputM);
y=y.*conj(y)/N;
subplot (3,2,2); plot (k,y); 
title ('FrFFT = FFT'); xlabel ('Frequency');

% alpha = 0.5 fractional FFT
alpha = 1/(2*N);
W=exp(-i*2*pi*alpha);

k=0:(M-1);
z=A.*W.^(-k);
subplot (3,2,3); plot (z,'*'); 
xlim([-1 1]); ylim([-1 1]);
xlabel('Real(z_n)'); ylabel('Imaginary(z_n)')
title ('\alpha = 1/2N, N=M=2^5')

display ('Fractional FFT: alpha = 1/2N')
tic; y = frfft ( X, alpha); toc;
outputM = length(y) 
k=linspace(0,N^2*alpha, outputM);
y=y.*conj(y)/N;
subplot (3,2,4); plot (k,y); 
title ('FrFFT'); xlabel ('Frequency');

% alpha = 1/3 fractional FFT
alpha = 1/(3*N);
W=exp(-i*2*pi*alpha);

k=0:(M-1);
z=A.*W.^(-k);
subplot (3,2,5); plot (z,'*'); 
xlim([-1 1]); ylim([-1 1]);
xlabel('Real(z_n)'); ylabel('Imaginary(z_n)')
title ('\alpha = 1/3N, N=M=2^5')

display ('Fractional FFT: alpha = 1/3N')
tic; y = frfft ( X, alpha); toc;
outputM = length(y) 
k=linspace(0,N^2*alpha, outputM);
y=y.*conj(y)/N;
subplot (3,2,6); plot (k,y); 
title ('FrFFT'); xlabel ('Frequency');



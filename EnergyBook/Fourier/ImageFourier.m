function [] = ImageFourier( N )
%ImageFourier plots the N harmonics of the complex
%exponential as well as the corresponding real cosine
%and imaginary sine components
%ImageFourier plots the complex exp(-i*2*pi/N) 
%for N = 1,...,16
%

close all
clc
set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold') 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

if (nargin ==0)
    N=8;
end
w=exp(-i*2*pi/N)
m=0:N-1;
mm=m'*m;
W=w.^(m'*m)%+1e-15*i
W(:,1)=W(:,1)+1e-15*i
for k=1:N
    subplot (3,N,k); plot (W(:,k),'-o')
    ylim([-1 1]);  xlim([-1 1]);
    if (k==1)
    ylabel('yi');
    end
    xlabel ('x');
    kstr=num2str(k-1);
    Nstr=num2str(N);
    tstr=['exp(-i2\pi' kstr '/' Nstr ')'];
    title(tstr);

    subplot (3,N,k+N);
    t=0:.01:N-1;
    plot (m,real(W(:,k)),'*',t,cos(2*pi*(k-1)/N.*t)); %hold on;
    tstr=['cos(n2\pi' kstr '/' Nstr ')'];
    title(tstr); xlabel ('n'); axis tight
    ylim([-1 1])
    if (k==1)
    ylabel('x');
    end

    subplot (3,N,k+2*N);
    plot (m,imag(W(:,k)),'*',t,-sin(2*pi*(k-1)/N.*t)); %hold on;
    tstr=['(-i)sin(n2\pi' kstr '/' Nstr ')'];
    title(tstr);xlabel ('n'); axis tight
    ylim([-1 1])
    if (k==1)
    ylabel('yi');
    end   
end

figure
for N=1:16
    w=exp(-i*2*pi/N);
    m=0:N-1;
    mm=m'*m;
    W=w.^(m'*m);
    if (N==1)
        W=W+1e-15*i;
    end
    %W=fft(eye(N));
    subplot (4,4,N);plot (W,'-o'); 
    ylim([-1 1]); xlim([-1 1])
    Nstr=num2str(N);Tstr=['exp(-i2\pikn/' Nstr ')'];
    title (Tstr);
end
end

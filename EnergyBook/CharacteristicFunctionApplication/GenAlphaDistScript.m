%GenAlphaDistScript creates several plots with random data
%vectors created by function GenAlphaDist
close all; clear all; clc
set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold') 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

sigma=0.01; mu=0; beta=0; 
loop=3

%Compare PDF for alpha=1, 1.5, 2
alphai=linspace(1,2,loop)
N=2^16 
for i=1:loop;
    alpha=alphai(i);
    Y= GenAlphaDist (alpha, beta, mu, sigma,  N );
    [yout,xout] = hist(Y,linspace(-10*sigma,10*sigma,200));
    y(i,:)=yout/N;
end
figure
subplot (2,1,1);
plot (xout,y(1,:),xout,y(2,:),'--',xout,y(3,:),':');
legend('\alpha=1','\alpha=1.5','\alpha=2')
ylabel('PDF(X)'); xlabel('X');
subplot (2,1,2); 
semilogy (xout,y(1,:),xout,y(2,:),'--',xout,y(3,:),':');
legend('\alpha=1','\alpha=1.5','\alpha=2')
ylabel('PDF(X)'); xlabel('X');

%Compare PDF for beta=-1, 0, 1
alpha=1.5
betai=linspace(-1,1,loop)
N=2^16 
for i=1:loop;
    beta=betai(i);
    Y= GenAlphaDist (alpha, beta, mu, sigma,  N );
    [yout,xout] = hist(Y,linspace(-10*sigma,10*sigma,200));
    y(i,:)=yout/N;
end
figure
subplot (2,1,1);
plot (xout,y(1,:),xout,y(2,:),'--',xout,y(3,:),':');
legend('\beta=-1','\beta=0','\beta=1')
ylabel('PDF(X)'); xlabel('X');
subplot (2,1,2); 
semilogy (xout,y(1,:),xout,y(2,:),'--',xout,y(3,:),':');
legend('\beta=-1','\beta=0','\beta=1')
ylabel('PDF(X)'); xlabel('X');

%Compare Simulated Return over time for alpha=1, 1.5, 2
beta=0;
alphai=linspace(1,2,loop)
N=2^8
for i=1:loop;
    alpha=alphai(i);
    Y= GenAlphaDist (alpha, beta, mu, sigma,  N );
    S(i,:)=100*exp(cumsum(Y));
end
figure
plot (1:N,S(1,:),1:N,S(2,:),'--',1:N,S(3,:),':');
ylabel('Price')
xlabel('Time ')
legend('\alpha=1','\alpha=1.5','\alpha=2')
axis tight



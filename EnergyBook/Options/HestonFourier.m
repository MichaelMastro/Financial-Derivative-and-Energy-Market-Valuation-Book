function [ c0 ] = HestonFourier (S0,K0,r,t,T,...
    kappa,lambda,sigma,rho,eta,v0)
% HestonFourier performs Fourier inversion of Heston model
% The option (or Greek) value is decoupled into a volatility 
% independent payoff function in Fourier Space
% (e.g., European Call: W(tau=0)=K.^(1+i*w)./(i*w-w.^2))
% and a Green Function, G=exp(C+vol*D), which is Fundamental 
% Transform if G(tau=0)=1 and G(w,vol,tau) satisfies 
% Heston style PDE. 
% Requires Numerical Integration of Solution
% C=exp(-r*tau)/2pi * integral

% HestonFourier will generate various plots comparing
% call price, Greeks, and implied volatility 
% vs. time to expiration, strike price, initial volatility
% and asset price

close all;
clc;
set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold'); 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

disp ('Heston Stochastic Volatility Option Calculation')
disp ('via Fourier inversion approach')

if (nargin < 10), v0 = 0.01; end
if (nargin < 9), eta = 0.3; end
if (nargin < 8), rho = -0.8; end
% usually stock price and volatility negative correlation,
% which leads to BS implied volatility smirk
if (nargin < 7), sigma = 0.5; end
if (nargin < 6), lambda =0.05; end 
if (nargin < 5), kappa = 1; end
if (nargin < 4), T =1; end
if (nargin < 4), t =0; end
if (nargin < 3), r = 0.03; end
if (nargin < 2), K0 = 50; end
if (nargin < 1), S0 = 50; end

tau=T-t;
prefactor=(1/(2*pi))*exp(-r*tau);
range=200 %+/- integration limits
% number of data points: strike, stock price and/or time
MaxCount=10 %%%%%%%%%%%%%%%%%%%%
OffSetImag=1.5; % offset contour of integration 
% call option: Fourier Transform Payoff=K0
jmid=round(MaxCount/2);

str1=['\kappa=' num2str(kappa)];
str2=[' \sigma=' num2str(sigma)];
str3=[' \lambda=' num2str(lambda)];
str4=[' \eta=' num2str(eta)];
str5=[' \rho=' num2str(rho)];
str6=[' v_0=' num2str(v0)];
textstr = [str1 str2 str3 str4 str5];
textstrplusv0 = [str1 str2 str3 str4 str5 str6];

disp ('Price vs. Asset Price and vs. Volatility');
disp('j k S(j,k) tau(j,k)');
for j=1:MaxCount
    for k=1:MaxCount
        S(j,k)=(2*K0/MaxCount)*j; 
        v(j,k)=0.2*k;
        disp([j k S(j,k) v(j,k)]);
        x=log(S(j,k))+r*tau;
               
        cInt=@(w) (exp(-i*w*x)...
            .*(K0.^(1+i*w)./(i*w-w.^2))...
            .*G(w,kappa,lambda,sigma, rho, eta, tau,v(j,k)));
        c(j,k)=real(prefactor*...
            quad(cInt,-range+2i,range+2i));
        
        deltaInt=@(w) ((-i.*w./S(j,k)).*exp(-i*w*x)...
            .*(K0.^(1+i*w)./(i*w-w.^2))...
            .*G(w,kappa,lambda,sigma, rho, eta, tau,v(j,k)));
        Delta(j,k)=real(prefactor*...
            quad(deltaInt,-range+2i,range+2i));
        
       % gammaInt=@(w) ((-w.^2/S(j,k)^2).*exp(-i*w*x)...
       %     .*(K0.^(1+i*w)./(i*w-w.^2))...
       %     .*G(w,kappa,lambda,sigma, rho, eta, tau,v(j,k)));
       % Gamma(j,k)=real(prefactor*quad(gammaInt,-range+2i,range+2i));
        
       % vegaInt=@(w) ((-w.^2/S(j,k)^2).*exp(-i*w*x)...
       %     .*(K0.^(1+i*w)./(i*w-w.^2))...
       %     .*VV(w,kappa,lambda,sigma, rho, eta, tau,v(j,k)));
       % Vega(j,k)=real(prefactor*quad(vegaInt,-range+2i,range+2i));
    end   
end

figure; subplot (2,1,1); surf(S,v,c); 
title(['Heston Inversion: Strike = ' num2str(K0) ', T = ' num2str(T)])
zlabel ('Call Value'); xlabel ('Asset Price'); 
ylabel ('Variance'); axis tight

subplot (2,1,2); surf(S,v,Delta); 
title(textstr); axis tight
zlabel ('Delta'); xlabel ('Asset Price'); ylabel ('Variance');

%figure; subplot (2,1,1); surf(S,v,Gamma); 
%title(['Fourier Inversion of Heston Model: K = ' num2str(K0)])
%zlabel ('Gamma'); xlabel ('Asset Price'); ylabel ('Variance');
%axis tight

%subplot (2,1,2);surf(S,v,Vega); title(textstr)
%zlabel ('Vega'); xlabel ('Asset Price'); ylabel ('Variance');
%axis tight
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c0=c(jmid,jmid); % call price returned from function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% call option value vs. initial asset price and vs. time
clear tau
tauk=linspace(0.1,5,MaxCount);
disp ('Price vs. Initial Asset Price and vs. Time');
disp('j k S(j,k) tau(j,k)');
for j=1:MaxCount
    for k=1:MaxCount       
        S(j,k)=(2*K0/MaxCount)*j;
        tau(j,k)=tauk(k);
        x=log(S(j,k))+r*tau(j,k); 
        disp([j k S(j,k) tau(j,k)]);
        cInt=@(w) (exp(-i.*w.*x)...
            .*(K0.^(1+i.*w)./(i.*w-w.^2))...
            .*G(w,kappa,lambda,sigma, rho, eta, tau(j,k),v0));
        prefactor=(1/(2*pi))*exp(-r*tau(j,k));
        c(j,k)=real(prefactor*...
            quadl(cInt,-range+2i,range+2i));
    end
end
figure; surf(S,tau,c ); 
title(['Heston: Strike K = ' num2str(K0) ' ' textstrplusv0])
zlabel ('Call Value'); xlabel ('Asset Price'); 
ylabel ('T_{expiration}'); axis tight
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate and plot Heston call option value and 
% Black-Scholes Implied volatility 
% vs. strike price and vs. time
% Negative correlation (rho) btw stock price and volatility in 
% Heston model generates a skew in the implied volatility
clear tau
tauk=linspace(0.1,5,MaxCount);
disp ('Calculate Call Price vs. Strike and vs. Time');
disp('j k K(j,k) tau(j,k)');
for j=1:MaxCount
    for k=1:MaxCount       
        K(j,k)=(2*K0/MaxCount)*j;
        tau(j,k)=tauk(k);
        x=log(S0)+r*tau(j,k);        
        
        disp([j k K(j,k) tau(j,k)]);
        cInt=@(w) (exp(-i.*w.*x)...
            .*(K(j,k).^(1+i.*w)./(i.*w-w.^2))...
            .*G(w,kappa,lambda,sigma, rho, eta, tau(j,k),v0));
        prefactor=(1/(2*pi))*exp(-r*tau(j,k));
        c(j,k)=real(prefactor*...
            quadl(cInt,-range+2i,range+2i));
    end
end
figure
subplot (2,2,1); surf(K,tau,c ); 
title(['Heston: S_0 = ' num2str(S0)])
zlabel ('Call Value'); xlabel ('Strike Price'); 
ylabel ('T_{expiration}');axis tight 

disp ('Calculate Implied Volatility');
% Find implied volatility with Matlab's fminsearch via function
% BSdiff, which compares Black-Scholes price 
% (at BS implied volaility) to Heston call price
volguess=v0+sigma;
d=0;
for j=1:MaxCount
  for k=1:MaxCount
    % fixed parameters 
    param = [c(j,k) K(j,k) S0 tau(j,k) r d]; 
    ImpVol(j,k) =...
        fminsearch(@(vol) BSdiff(vol,param),volguess);
    volguess=ImpVol(j,k); % speeds up next fminsearch
  end
end
subplot (2,2,3); surf(K,tau,ImpVol ); title(textstrplusv0);
zlabel ('Implied Volatility'); xlabel ('Strike Price'); 
ylabel ('T_{expiration}'); axis tight

% Take a slice at first time point 
subplot (2,2,2); plot(K(:,1), c(:,1) )
ylabel ('Call Value'); xlabel ('Strike Price'); 
title(['T = ' num2str(tau(1,1)) ' Years']);
axis tight

subplot (2,2,4); plot(K(:,1), ImpVol(:,1) )
ylabel ('Implied Volatility'); xlabel ('Strike Price'); 
title(['T = ' num2str(tau(1,1)) ' Years']);
axis tight

end  %%%%%%%%%%%%%%%% End HestonFourier Function%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Gval D] = G(w,kappa,lambda,sigma, rho, eta , tau,v)
% G is Integrand for call price, Delta, Gamma equations
d=sqrt( (w.^2-i*w)*(sigma^2) + (kappa+lambda+i*sigma*rho*w).^2);
g=(kappa+lambda+ i.*sigma.*rho.*w +d)./...
    (kappa+lambda+ i.*sigma.*rho.*w -d);
D=( (kappa+lambda+ i.*sigma.*rho.*w+d).*(1-exp(d.*tau)) ) ./...
       ( sigma.^2 .*  (1-g.*exp(d.*tau)) );

 %   C=(kappa.*eta./sigma.^2)...
 %   .* ( (kappa+lambda+i.*sigma.*w + d).*tau...
 %   -2.*log( (1-g.*exp(d.*tau))/(1-g)));

 % Analytical Equation for C can be unstable, therefore use
 % Brute Force Integration -> C=kappa*eta* integral(D dt)
   maxcounter=50;
   Dsum=0;
for l=1:maxcounter
   counter=(l-1)/(maxcounter-1);
    Dsum= Dsum+ (( (kappa+lambda+ i.*sigma.*rho.*w+d)...
        .*(1-exp(d.*tau.*counter)) ) ./...
       ( sigma.^2 .*  (1-g.*exp(d.*tau.*counter)) ));
end
C=kappa*eta*Dsum/(maxcounter-1);
Gval=exp(C + D*v);
end

function VegaVal = VV(w,kappa,lambda,sigma, rho, eta , tau,v)
% VegaVal Provides the Integrand for integration in 
% equation for Vega. Expression has extra D in 
% VegaVal=D.*exp(C + D*v); as compared to equations
% for call price, Delta, Gamma
% Easier to use separate function
[Gval D] = G(w,kappa,lambda,sigma, rho, eta , tau,v);
VegaVal=D.*Gval;
end

function diff2 = BSdiff(vol,param)
% BSdiff returns squared difference of Black-Schole price
% for input volatility to input call price (from Heston
% model in this application)
c=param(1);
K=param(2);
S=param(3);
T=param(4);
r=param(5);
d=param(6);

diff=c-BlackScholesCall (K,S,T,vol,r,d);
diff2=diff^2;
end

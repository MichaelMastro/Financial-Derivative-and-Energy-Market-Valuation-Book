function C = WhaleyFuture ( );
%WhaleyFuture Calculates Call Value of Option on a Futures 
%Derived for an asset with cost of carry b in 
%Barone-Adesi and Whaley [1987] J. of Finance.
%Futures Contract has a cost of carry of zero, 
%thus equations can be simplified by setting b=0
close all;
clc;
set(0,'defaultaxeslinewidth',3); set(0,'defaultlinelinewidth',3);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold'); 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);


tic
r=0.3;vol=0.3;K=50;

deviation=1000;
Fstar=K*1.5;
T=1;
%Newton-Raphson Search to find critical Futures Price (Fstar)
%where indifferent to exercising (F-K) vs. holding American call
%Newton-Raphson has quadratic approach thus very small deviation
%adds minimal additional iterations

Mdivh=2*r/(vol^2.*(1-exp(-r*T))); % Pre-calculate since not
q2= 0.5*(1+sqrt(1+4.*Mdivh));     % dependent on Fstar
while (deviation>1e-10)
    d1 = ( log(Fstar./K)+ (vol^2/2).*T)./ (vol.*sqrt(T));
    g=BlackCall (K,Fstar,T,vol,r)...
        +(Fstar/q2).*(1-exp(-r*T)*myNormCDF(d1))...
        -Fstar+K;
    gprime=-1+(exp(-r*T)*myNormCDF(d1))*(1-1/q2)...
        +(1/q2)*(1-(exp(-r*T)*myNormPDF(d1))/(vol*sqrt(T)));
    deviation=g/gprime;
    Fstar=Fstar-deviation;
    deviation=abs(deviation);
end
toc
fprintf(1, 'Critical Futures Price = %6.2f \n', Fstar); 

Fmin=1; Fmax=100; Fstep=0.1;
FstarIndex=round((Fstar-Fmin)/Fstep)+1;
F=Fmin:Fstep:Fmax; 

cEuropean=BlackCall (K,F,T,vol,r);
%Vectorize Code and avoid for-if loop to speed up code
CAmerican=zeros(1,length(F));
%C=Quadrature Approx  F<Fstar
A2(1:FstarIndex-1)=(Fstar/q2)*(1-exp(-r*T)*myNormCDF(d1));
CAmerican(1:FstarIndex-1)=cEuropean(1:FstarIndex-1)...
    +A2(1:FstarIndex-1).*(F(1:FstarIndex-1)./Fstar).^q2;
%C=F-K  F>=Fstar
CAmerican(FstarIndex:end)=F(FstarIndex:end)-K;

FK=max(0,F-K);
FKemrt=exp(-r*T)*FK;
FstarK=Fstar-K;

plot (F,CAmerican,F,cEuropean,'-.',F,FK,'--',F,FKemrt,':',...
    Fstar,FstarK,'+',K,0,'*');
xlabel('Futures Price'); ylabel('Call Value');
legend ('American call','European call','F-K','e^{-rT}(F-K)',...
    'Location','NorthWest');
text(K-1,-1,'K');
text(Fstar-2,FstarK+2,'F*');
text(K+15,10,'e^{-rT}(F-K)')
text(K+25,40,'(F-K)')
title('Barone-Adesi and Whaley Quadratic Approximation');   

end

function c = BlackCall (K,F,T,vol,r);
    d1 = ( log(F./K)+ (vol^2/2).*T)./ (vol.*sqrt(T));
    d2 = d1-vol*sqrt(T);
    c = exp(-r*T).*(F.*myNormCDF(d1)-K.*myNormCDF(d2));
end

function p = BlackPut (K,F,T,vol,r);
    d1 = ( log(F./K)+ (vol^2/2).*T)./ (vol.*sqrt(T));
    d2 = d1-vol*sqrt(T);
    p = exp(-r*T).*(K.*myNormCDF(-d2)-F.*myNormCDF(-d1));
end


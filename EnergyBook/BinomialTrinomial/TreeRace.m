function TreeRace()
%TreeRace compares Trinomial v. Binomial Tree calcs.
close all
clc

set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold') 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

S0=50;
K=50;
r=0.1;
sigma=0.5;
T=1;
d=0; %dividend yield
CallFlag=0;
EuropeanFlag=0; 

trials=16;
TriValues=zeros(1,trials);
TriTime=zeros(1,trials);

for n=1:trials;
    Ntrial(n)=(n+2)^2;
    tic
    TriValues(n)=VectorTrinomial(S0,K,r,sigma,T,Ntrial(n),d,...
                        CallFlag,EuropeanFlag);
    TriTime(n)=toc;
end

for n=1:trials;
    tic
    BiValues(n)=VectorBinomial(S0,K,r,sigma,T,Ntrial(n),d,...
                        CallFlag,EuropeanFlag);
    BiTime(n)=toc;
end

figure
subplot(3,1,1)
semilogy(Ntrial,abs(TriValues-BiValues))
title('Binomial vs. Trinomial Tree Comparison'); 
ylabel('log(\DeltaV)');
axis tight

subplot(3,1,2)
plot(Ntrial,TriValues,Ntrial,BiValues)
ylabel('Value');
if ((CallFlag==1) && (EuropeanFlag==1))
    title('European Call'); 
elseif ((CallFlag==1) && (EuropeanFlag==0))
    title('American Call'); 
elseif ((CallFlag==0) && (EuropeanFlag==0))
    title('American Put');
else
    title('European Put');
end  
axis tight

subplot(3,1,3)
plot(Ntrial,TriTime,Ntrial,BiTime)
legend('Trinomial','Binomial','location','NorthWest');
xlabel('Steps');ylabel('Calculation Time');
                    
end

% GammaRandScript with GammaRand simulate Gamma process
% Creates one plot with four subfigures of
% delta movement vs time
% total (summed) movement vs time
% linear PDF of movement
% log PDF of movement

close all; clear all; clc
set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold') 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

a=10;  
b=10; 
dt=0.01 
adt=a*dt; %a*dt=k
t=0:dt:1;
lengtht=length(t);

for j=1:lengtht
	Gab(j)=GammaRand(adt,b);
end
SumGab=cumsum(Gab);
disp('min(Gab)');disp(min(Gab));
disp('max(Gab)');disp(max(Gab));
GjumpSize=linspace(min(Gab), max(Gab), lengtht);

pdfGamma=b.^adt.* GjumpSize.^(adt-1).*exp(-GjumpSize*b)...
            ./ (  gamma(adt) ); 
% normalize continuous function
pdfGamma=pdfGamma/sum(pdfGamma(3:end)); 

figure
subplot (4,1,1)
plot (t,Gab)
title('Gamma Process')
xlabel('Time')
ylabel('\Delta')

subplot (4,1,2)
plot (t, SumGab)
xlabel('Time')
ylabel('\Sigma\Delta')

subplot (4,1,3)
[n,xout]=hist (Gab,round(lengtht/4)); %
bar (xout, n/sum(n))
hold on
plot (GjumpSize, pdfGamma,'g--') 
ylim([0, 1]) 
xlabel('\Delta')
ylabel('PDF(\Delta)')

subplot (4,1,4)
semilogy (xout, n/lengtht)
hold on
semilogy (GjumpSize, pdfGamma,'g--') 
ylim([0.001, 1]) 
xlabel('\Delta')
ylabel('log(PDF(\Delta))')
legend('Discrete', 'Continuous')

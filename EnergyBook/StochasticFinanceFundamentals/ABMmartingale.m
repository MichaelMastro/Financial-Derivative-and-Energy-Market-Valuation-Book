%ABMmartingale models stock price movement of ABM:
%St-S0 ~ normal(mu*t, sigma*t^0.5)
%xi function tranform probability density to create mertingale
close all
clear all
clc
set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold');    
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

mu=0.8;
t=1;
sigma=0.4;
delS=[-4:0.1:4];

pdfP = myPDF (delS, mu*t, sigma*sqrt(t));
pdfPmod = myPDF (delS, 0, sigma*sqrt(t));
XI=exp((-1/sigma^2)*(mu*delS-0.5*mu^2*t));
semilogy (delS, pdfP, '--',  delS, XI, delS, pdfPmod,'-.')
legend ('PDF' ,'\xi', 'modified PDF', 'Location', 'South')
xlabel('\Delta S_t'); ylabel('Probability Density')%%
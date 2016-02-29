% MertonScript shows increase in call value of Merton
% with jumps vs. Black-Scholes without jumps
% underlying Black-Scholes volatility assumed equal

close all; clear all; clc
set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold') 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

T=1;        Ttext=['T = ' num2str(T)];    
r=0.05;     rtext=['r = ' num2str(r)];    
K=100;      Ktext=['K = ' num2str(K)];    
%S0=100;     Stext=['S_0 = ' num2str(S0)];    

%%% Black-Scholes Parameters %%%
volBS = 0.2; Voltext=['\sigma_{BS} = ' num2str(volBS)];    
d=0;

%%% Merton Series %%%
muJ = 0.1;       
sigmaJ = 0.1;    
lambda = 0.5;    
mmuJtext=['\mu_{Jump} = ' num2str(muJ)];    
msigmaJtext=['\sigma_{Jump} = ' num2str(sigmaJ)];    
lambdatext=['\lambda = ' num2str(lambda)];    

Srange=1:180;
for S=Srange
    cBS(S) = BlackScholesCall (K,S,T,volBS,r,d);
    cM(S) = MertonSeriesCall (K,S,T,volBS,...
                        r, d, muJ, sigmaJ, lambda);    
end
plot (Srange, cM, Srange, cBS,'--')
axis tight
legend ('Merton','Black-Scholes','Location','NorthWest')
ylabel('Call Price')
xlabel('S_0')

bsStr(1)= {'Black-Scholes'};
bsStr(2)= {rtext};
bsStr(3)= {Voltext};
text(125, 20, bsStr)

genStr(1)= {Ktext};
genStr(2)= {Ttext}; 
% genStr(3)= {Stext};
text(90, 60, genStr);

JumpStr(1)= {'Merton '};
JumpStr(2)= {mmuJtext};
JumpStr(3)= {msigmaJtext};
JumpStr(4)= {lambdatext};
text(50, 20, JumpStr);

title ('Merton Log-Normal Jump Process vs. Black-Scholes')

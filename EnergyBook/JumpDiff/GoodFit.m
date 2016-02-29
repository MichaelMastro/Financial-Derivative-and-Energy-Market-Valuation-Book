function GoodFit(f1,f2,n, binCenters)
% GootFit function provides Goodness of Fit measure
% assume two CFD vectors are the same length
CDF1=cumsum(f1)/n;
CDF2=cumsum(f2)/n;
% Kolmogorov–Smirnov statistic 
Dn1=max(abs(CDF1-CDF2));
dn1=Dn1*sqrt(n);
% dn=scaled distances between two CDFs
j=1:1000; %a ssume 1000 ~ infinity for this calculation
Pd1=signif(dn1); % Pd1=significance level of scales value of dn
%
fprintf(1,'Kolmogorov–Smirnov \n');
fprintf(1,...
'Dn = %6.4f Scaled dn = %6.4f Significance = %6.4f \n',...
Dn1, dn1, Pd1);

d=0.01:0.01:2.5;
pd= signif(d);

d05band=1.36/sqrt(n); 
fprintf(1,'D-alpha = %6.4f \n',d05band);

d05up=CDF1+d05band;
d05down=CDF1-d05band; 
figure
subplot (1,2,1); plot (binCenters, d05up, '-.' ,...
    binCenters, CDF1, binCenters,...
    CDF2, binCenters, d05down, '-.' )
legend('+d_{\alpha=0.05}','\Phi_{JD}','\Phi_{Data}',...
    '-d_{\alpha=0.05}','location','NorthWest')
axis tight
xlim([-0.05, 0.05])
title('Reject CDF if Empirical CDF falls outside +/-d_{\alpha=0.05}')
text(0 ,0.67,'\downarrow');
text(0 ,0.5,'\uparrow 5% Critcal Level');
xlabel('Log-Return'); ylabel('Cumulative Distribution')

subplot (1,2,2); plot (d,pd)
xlabel('Scaled Max CDF Error d=n^{0.5}D_n'); 
ylabel('Significance \alpha ')
text(dn1 ,Pd1,'\downarrow JD Significance Level',...
    'VerticalAlignment','bottom')
title('Kolmogorov–Smirnov Significance Curve')

% H0 CFD1=CFD2 confidence 
% (no basis to reject at significance level alpha)
% if dn<q(1-alpha)

% Cramér–von-Mises 
% Wn2=n*sum((CFD1-CFD2).^2);
% fprintf(1, 'Cramér–von-Mises Wn2 = %6.4f \n',  Wn2);
end

function Pd= signif (d);
ld=length(d);
j=1:1000; % assume 1000 ~ infinity for this calculation
for i=1:ld
    Pd(i)=2*sum((-1).^(j-1).*exp(-2.*j.^2.*d(i).^2));
end
end
function ProbPlot ()
%ProbPlot examines the probabilities in Hull-White 
%Mean Reverting Tree. Hull-White select the maximum tree
%size and thus the location of the edge to avoid
%a negative pm in the interior at abs(jM)>0.815
%But a negative probability at the top or bottom
%edge can exist if abs(jM)<0.185

close all
clc

set(0,'defaultaxeslinewidth',3); set(0,'defaultlinelinewidth',3);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold') 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

sigma=0.2;
Steps=10;
T=1;
a=0.1;

dt=T/Steps
M=-a*dt %=exp(-a*dt)-1;
V=sigma^2*dt;%=sigma^2*(1-exp(-2*a*dt))/(2*a);

nodes=10*ceil(-0.1835/M);
jtotal=2*nodes+1

jmid=nodes+1

nodes=nodes;
jrange=(jmid-nodes):1:(jmid+nodes);

for j=jrange    
    jreal=jmid-j; 
    %Top
    if (jreal>0)
        puTop(j)=7/6+(jreal^2*M^2+3*jreal*M)/2;
        pmTop(j)=-1/3-jreal^2*M^2-2*jreal*M;
        pdTop(j)=1/6+(jreal^2*M^2+jreal*M)/2;
    else
       puTop(j)=NaN; pmTop(j)=NaN;pdTop(j)=NaN;
    end
    
    %bottom
    if (jreal<0)
    puBottom(j)=1/6+(jreal^2*M^2-jreal*M)/2;
    pmBottom(j)=-1/3-jreal^2*M^2+2*jreal*M;
    pdBottom(j)=7/6+(jreal^2*M^2-3*jreal*M)/2;
    else
       puBottom(j)=NaN; pmBottom(j)=NaN;pdBottom(j)=NaN;
    end
        
    %Interior
    puMid(j)=1/6+(jreal^2*M^2+jreal*M)/2;
    pmMid(j)=2/3-jreal^2*M^2;
    pdMid(j)=1/6+(jreal^2*M^2-jreal*M)/2;

end% j loop
    
jmax=ceil(-0.1835/M)*-1*M
jR=(jmid-jrange)*-1*M;
PlusOne=ones(length(jR),1);
MinZero=zeros(length(jR),1);

figure
subplot (3,1,1)
plot(jR,puTop,':',jR,pmTop,jR,pdTop,'--',...
    jR,PlusOne,'-',jR,MinZero);
hold on
stem(jmax,1)
title('Top Layer Nodes')
legend ('p_u','p_m','p_d','location','West')
axis tight
hold off

subplot (3,1,2)
plot(jR,puMid,':',jR,pmMid,jR,pdMid,'--',...
    jR,PlusOne,'-',jR,MinZero);
hold on
stem([-jmax jmax],[1 1]);
title('Interior Nodes')
ylabel('Transition Probability');
axis tight
hold off

subplot (3,1,3)
plot(jR,puBottom,':',jR,pmBottom,jR,pdBottom,'--',...
    jR,PlusOne,'-',jR,MinZero);
title('Bottom Layer Nodes')
legend ('p_u','p_m','p_d','location','East')
hold on
stem(-jmax,1);
axis tight
xlabel('-jM');
hold off

end

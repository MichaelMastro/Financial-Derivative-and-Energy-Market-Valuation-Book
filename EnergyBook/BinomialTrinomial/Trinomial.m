function Trinomial (ReadDataFlag)
%Trinomial Reads/Defines data and calls HullWhiteTrinomial
%if ReadDataFlag==1 then trinomial tree is fit to current
%futures curve (as of January 2011); and sigma and drift
%are extracted from historical WTI front futures contract
%if ReadDataFlag==0 (or anything else including []) then
%replicate Hull, Options, Futures and Other Derivatives
%6th ed, p 720 example
%Using mean reversion process:
    %dS/S=kappa[mu-lnS]dt+sigma*dz
    %dX=kappa[alpha-X]dt+sigma*dz
    %alpha=mu-sigma^2/(2*kappa)
%Estimate parameters from Historic Futures Price data
%via a Least Squares Fit
%An alternate but more involved process is to extract 
%market participants implied sigma, etc. from 
%current Futures Options prices

close all
clc
set(0,'defaultaxeslinewidth',3); set(0,'defaultlinelinewidth',3);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold') 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

if ((nargin == 1) && (ReadDataFlag ==1))
    %read historical front WTI future contract 
    %1983 to 2011-1-11
    HistoricalPrice=load('WTIfrontContract.dat'); 
    %Input Matrix should have 2 columns: time , price
    THist=(HistoricalPrice(:,1));
    FHist=(HistoricalPrice(:,2));
    figure; plot (1983+4/12+THist,FHist);
    
    lnFHist=log(FHist);
    HistTimeStep=HistoricalPrice(2,1)-HistoricalPrice(1,1);
    %dS/S=kappa(mu-lnS)dt+sigma*dz
    %LeastSquaresOU = Mean Reversion %using X=lnS 
    %dX=kappa(alpha-X)dt+sigma*dz
    %alpha=mu-sigma^2/(2*kappa)
    [alpha,kappa, sigma] = LeastSquaresOU (lnFHist,HistTimeStep)

    %in the structure of Hull-White 
    %theta=kappa*mu and a=kappa
    %d(lnS)=[theta-a*lnS]dt*sigma*dz
    %To construct the 1st step driftless tree, temporarily 
    %assume theta=0, X=lnS
    %dX=-a*dt+sigma*dz
    r=0.05
    Steps=50;
    T=3;
    a=kappa;    %and sigma=sigma
    
    %read directly crude oil wti future data 2011-1-11
    FuturesPrice=load('CrudeOilFuturePrice.dat'); 
    %Input Matrix should have 2 columns: time , price
    Tvector=(FuturesPrice(:,1));
    Fvector=(FuturesPrice(:,2));
    xlabel('Time [Years]'); ylabel('Crude Oil Futures Price');
    title('Historical WTI Front Futures Contract'); 
    axis tight
else
    %Replicate Hull example in Options, Futures and Other...
    %Derivatives, 6th Ed.,p720
    Tvector=[0,1,2,3]
    Fvector=[20,22,23,24]
    sigma=0.2;
    Steps=3;
    T=3;
    a=0.1;
    kappa=a;
    r=0.05;
end
    
[S, pu, pm, pd, Arrow, jmax, kmatrix] = HullWhiteTrinomial...
    (Tvector,Fvector,T,Steps,sigma,a, r);

K=Fvector(1)
tOption=2.6
ForwardT=3.5

fprintf(1,'\nEuropean Option Pricing via Trinomial Tree\n');
[EurCall, EurPut, EurCallFuture, EurPutFuture ]...
    =EurHWTrinomial (S, Arrow, T, Steps,...
    jmax, K, tOption,ForwardT,Tvector,Fvector,kappa, sigma)

fprintf(1,'\nAnalytical European Options via Blacks model\n');
[BlackCall, BlackPut ]...
    = BlackOptionMeanRev (K, tOption, tOption, ...
        r,Tvector,Fvector,kappa, sigma)

[BlackCallFuture, BlackPutFuture ]...
    = BlackOptionMeanRev (K, tOption, ForwardT, ...
        r,Tvector,Fvector,kappa, sigma)

fprintf(1,'\nAmerican Option Pricing via Trinomial Tree\n');
[AmerCall, AmerPut, AmerCallFuture, AmerPutFuture ]...
    =AmerHWTrinomial (S, pu, pm, pd, T, Steps,...
    jmax, K, tOption, kmatrix, r,...
    ForwardT,Tvector,Fvector,kappa, sigma) 

    
end

function [cBlack, pBlack ]...
    = BlackOptionMeanRev (K, tOption, ForwardT, ...
        r,Tvector,Fvector,kappa, sigma)
%BlackOptionMeanRev provides an analytical European Call and Put
%Valuation under mean reversion
%A Typical Futures Contract will have an option expiration (tOption)
%that precedes the Futures Contract maturity (ForwardT)
%Setting tOption=ForwardT will solve a spot option
%with expiration (tOption)

s=sqrt(sigma^2/(2*kappa)*exp(-2*kappa*(ForwardT-tOption))...
        *(1-exp(-2*kappa*tOption)));
FT = interp1(Tvector,Fvector,ForwardT, 'linear','extrap');
d1=log(FT/K)/s+0.5*s;
d2=d1-s;
%P=exp(-r*tOption) assumes constant interest rates
cBlack=exp(-r*tOption)*(FT*myNormCDF(d1)-K*myNormCDF(d2));
pBlack=exp(-r*tOption)*(K*myNormCDF(-d2)-FT*myNormCDF(-d1));


end

function [EurCall, EurPut, EurCallFuture, EurPutFuture ]...
    =EurHWTrinomial (S, Arrow, T, Steps,...
    jmax, K, tOption,ForwardT,Tvector,Fvector,kappa, sigma)
%EurHWTrinomial calculates European Call and Put
%at the expiration nodes using the Arrow-Debreu State prices

    dt=T/Steps;
    jtotal=2*jmax+1;

    IndexMaturity=round(tOption/dt)+1;
    %Assume large number of time steps thus a time level will be
    %very close to tOption. Better to increase time levels rather 
    %than use interpolation
    jmid=jmax+1; %disp(Arrow)
    nodes=min(IndexMaturity-1,jmax);
        jrange=(jmid-nodes):1:(jmid+nodes);
        EurCall=sum(Arrow(jrange,IndexMaturity)...
            .*max(S(jrange,IndexMaturity)-K,0) );
        EurPut=sum(Arrow(jrange,IndexMaturity)...
              .*max(K-S(jrange,IndexMaturity),0) );
 
%From Market Futures Curve Data Interpolate or Extrapolate 
%Future Price at time of Option and Time of Futures contract
%to be valued
    Ft = interp1(Tvector,Fvector,tOption)
    FT = interp1(Tvector,Fvector,ForwardT, 'linear','extrap')
 %Calculate Futures price at nodes using Clewlow and Strikland
 %formula in 'Valuing Energy Options in a One Factor...'
    FSt=FT*(S(jrange,IndexMaturity)/Ft)...
        .^exp(-kappa*(ForwardT-tOption)...
        *exp( (-sigma^2/(4*kappa)) * (exp(2*kappa*tOption)-1)...
          *(exp(-kappa*ForwardT)-exp(-kappa*tOption))));
    EurCallFuture=sum(Arrow(jrange,IndexMaturity)...
        .*max(FSt(jrange)-K,0) );
    EurPutFuture=sum(Arrow(jrange,IndexMaturity)...
          .*max(K-FSt(jrange),0) ); 
end

function [AmerCall, AmerPut,...
    AmerCallFuture, AmerPutFuture ]...
    =AmerHWTrinomial (S, pu, pm, pd, T, Steps,...
    jmax, K, tOption, kmatrix, r,...
    ForwardT,Tvector,Fvector,kappa, sigma) 
%AmerHWTrinomial calculates Put and Call by stepping backwards
%in time through a trinomial tree

    dt=T/Steps;
    jtotal=2*jmax+1;
    jmid=jmax+1;

    IndexMaturity=round(tOption/dt)+1;
    nodes=min(IndexMaturity-1,jmax);
    jrange=(jmid-nodes):1:(jmid+nodes);
    C(jrange,IndexMaturity)=max(S(jrange,IndexMaturity)-K,0);
    P(jrange,IndexMaturity)=max(K-S(jrange,IndexMaturity),0);

%Future Price at time of Option and Time of Futures contract
%to be valued
    Ft = interp1(Tvector,Fvector,tOption);
    FT = interp1(Tvector,Fvector,ForwardT, 'linear','extrap');

 %Calculate Futures price at nodes using Clewlow and Strikland
 %formula in 'Valuing Energy Options in a One Factor...'
    FSt(jrange,IndexMaturity)=FT*(S(jrange,IndexMaturity)/Ft)...
        .^exp(-kappa*(ForwardT-tOption)...
        *exp( (-sigma^2/(4*kappa)) * (exp(2*kappa*tOption)-1)...
          *(exp(-kappa*ForwardT)-exp(-kappa*tOption))));
    CFuture(jrange,IndexMaturity)=...
        max(FSt(jrange,IndexMaturity)-K,0);
    PFuture(jrange,IndexMaturity)=...
        max(K-FSt(jrange,IndexMaturity),0);

    emrdt=exp(-r*dt);
for i=IndexMaturity-1:-1:1
    nodes=min(i-1,jmax);
    jrange=(jmid-nodes):1:(jmid+nodes);
    
    FSt(jrange,i)=FT*(S(jrange,IndexMaturity)/Ft)...
        .^exp(-kappa*(ForwardT-tOption)...
        *exp( (-sigma^2/(4*kappa)) * (exp(2*kappa*tOption)-1)...
          *(exp(-kappa*ForwardT)-exp(-kappa*tOption))));
      
    for j=jrange
        k=kmatrix(j,i);
        C(j,i)=emrdt*(pu(j,i)*C(j-1+k,i+1)...
           +pm(j,i)*C(j+k,i+1)... 
            +pd(j,i)*C(j+1+k,i+1)); 
        C(j,i)=max(C(j,i),S(j,i)-K); %Amer Premium
        P(j,i)=emrdt*(pu(j,i)*P(j-1+k,i+1)...
           +pm(j,i)*P(j+k,i+1)... 
            +pd(j,i)*P(j+1+k,i+1)); 
        P(j,i)=max(P(j,i),K-S(j,i)); %Amer Premium
        
        CFuture(j,i)=emrdt*(pu(j,i)*CFuture(j-1+k,i+1)...
           +pm(j,i)*CFuture(j+k,i+1)... 
            +pd(j,i)*CFuture(j+1+k,i+1)); 
        %Amer Premium
        CFuture(j,i)=max(CFuture(j,i),FSt(j,i)-K); 
        
        PFuture(j,i)=emrdt*(pu(j,i)*PFuture(j-1+k,i+1)...
           +pm(j,i)*PFuture(j+k,i+1)... 
            +pd(j,i)*PFuture(j+1+k,i+1)); 
        %Amer Premium
        PFuture(j,i)=max(PFuture(j,i),K-FSt(j,i));        
    end
end
    AmerCall=C(jmid,1); AmerPut=P(jmid,1);
    AmerCallFuture=CFuture(jmid,1);
    AmerPutFuture=PFuture(jmid,1);
end
  

function [S, pu,pm, pd, Arrow, jmax, kmatrix]...
    = HullWhiteTrinomial(Tvector,Fvector,T,Steps,sigma,a, r)
%HullWhiteTrinomial fits a trinomial tree by Hull-White procedure
%Step 1: create trinomial tree of dX=-a*dt+sigma*dz
%where X=lnS 
%Step 2: fit nodes at each time step in trinomial tree
%to that futures price at that time step: F(t) = E[S(t)]
    dt=T/Steps
    M=-a*dt;%=exp(-a*dt)-1;
    V=sigma^2*dt;%=sigma^2*(1-exp(-2*a*dt))/(2*a);

    jmax=ceil(-0.1835/M)
    jtotal=2*jmax+1

    Q=zeros(jtotal,Steps+2);
    Arrow=zeros(jtotal,Steps+2);
    deltaX=sqrt(3*V);%sigma*sqrt(3*dt)
    %t=i*dt
    %X=j*deltaX;
    jmid=jmax+1
    X(jmid,1)=0;%=X0
    Q(jmid,1)=1;
    Arrow(jmid,1)=1;

for i=1:Steps+1
    %fprintf(1,'time step i=%6.0f\n',i);
    nodes=min(i-1,jmax);
    jrange=(jmid-nodes):1:(jmid+nodes);
    for j=jrange
        %fprintf(1,'j=%6.0f\t',j);
        if (j>1)
        X(j-1,i+1)=X(j,i)+deltaX;
        end
        X(j,i+1)=X(j,i);
        if (j<jtotal)
        X(j+1,i+1)=X(j,i)-1*deltaX;
        end
        
        jreal=jmid-j; 
        if (j==1) %top
            pu(j,i)=7/6+(jreal^2*M^2+3*jreal*M)/2;
            pm(j,i)=-1/3-jreal^2*M^2-2*jreal*M;
            pd(j,i)=1/6+(jreal^2*M^2+jreal*M)/2;
            k=1;  kmatrix(j,i)=k; %Simplify Amer Option
            
        elseif (j==jtotal) %bottom
            pu(j,i)=1/6+(jreal^2*M^2-jreal*M)/2;
            pm(j,i)=-1/3-jreal^2*M^2+2*jreal*M;
            pd(j,i)=7/6+(jreal^2*M^2-3*jreal*M)/2;
            k=-1; kmatrix(j,i)=k; %Simplify Amer Option
            
        else
            pu(j,i)=1/6+(jreal^2*M^2+jreal*M)/2;
            pm(j,i)=2/3-jreal^2*M^2;
            pd(j,i)=1/6+(jreal^2*M^2-jreal*M)/2;
            k=0; kmatrix(j,i)=k; %Simplify Amer Option           
        end

%Special case: Qs not discounted since const interest rates
        Q(j-1+k,i+1)=Q(j-1+k,i+1)+pu(j,i)*Q(j,i);
        Q(j+k,i+1)=Q(j+k,i+1)+pm(j,i)*Q(j,i);
        Q(j+1+k,i+1)=Q(j+1+k,i+1)+pd(j,i)*Q(j,i);  
    %Arrow-Debreu Prices discounted at the riskfree rate to the 
    %present value
        Arrow(j-1+k,i+1)=Arrow(j-1+k,i+1)...
            +pu(j,i)*Arrow(j,i)*exp(-r*dt);
        Arrow(j+k,i+1)=Arrow(j+k,i+1)...
            +pm(j,i)*Arrow(j,i)*exp(-r*dt);
        Arrow(j+1+k,i+1)=Arrow(j+1+k,i+1)...
            +pd(j,i)*Arrow(j,i)*exp(-r*dt);     
    end% j loop
    
    if (i>1)
        F(i) = interp1(Tvector,Fvector,(i-1)*dt);
        theta(i)=log(F(i)./ sum((Q(jrange,i)).*exp(X(jrange,i))));
    else
        F(1) = interp1(Tvector,Fvector,(i-1)*dt, 'linear','extrap');
        theta(1)=log(min(F(1),Fvector(1)));
    end
    S(jrange,i)=exp(theta(i)+ X(jrange,i)); 
end

    figure
    plot (Tvector,Fvector,0,F(1),'+','MarkerSize',10);
    xlabel('Time [Years]'); ylabel('Crude Oil Futures Price');
    title('WTI Futures as of January 2011');

    if (Steps<6)
        %disp('X');disp(X); %disp('Q');disp(Q);   
        disp('S');disp(S);%disp('theta');disp(theta);
        figure %Plot Tree for Spot Price of Oil
    for i=1:Steps
        nodes=min(i-1,jmax);
        %nodelength=nodes+1
        jrange=(jmid-nodes):1:(jmid+nodes);
        for j=jrange
           clear x y
            if (j==1) %top
                k=1;
            elseif (j==jtotal) %bottom
                k=-1;
            else   
                k=0;
            end
           x(:,1)=[(i-1)*dt i*dt];%+dt
           x(:,2)=x(:,1);  x(:,3)=x(:,1);
           y(:,1)=[S(j,i) S(j-1+k,i+1)]; %up
           y(:,2)=[S(j,i) S(j+k,i+1)]; %middle
           y(:,3)=[S(j,i) S(j+1+k,i+1)]; %down
           plot(x,y)
           txstr1(1) = {num2str(S(j,i),'S=%.1f')};
           txstr1(2) = {num2str(Q(j,i),'Q=%.2f')}; 
           text((i-1)*dt,S(j,i),txstr1);
           hold on
        end     
    end
    t=Steps+1;
    for j=1:jtotal  % Expiration nodes
       txstr2(1) = {num2str(S(j,t),'S=%.1f')};
       txstr2(2) = {num2str(Q(j,t),'Q=%.2f')}; 
       text((t-1)*dt*0.95,S(j,t),txstr2);
       hold on
    end
    %axis tight  
    xlabel('Time (Years)'); ylabel('Asset Value');
    hold off
    end %only plot if small number of time levels
end

function [alpha, sigma, kappa] = LeastSquaresOU (S,delta)
%LeastSquaresOU performs an  unweighted least squares
%fit to Orstein Uhlenbeck process
%Derivation in Press, Flannery, Teukolsky, Vetterling, 
%Numerical Recipes as well as an unweighted estimation  
%procedure by M.A. van den Berg available at www.sitmo.com

x=S(1:end-1);
y=S(2:end);

n= length (y);
xMean=mean(x);
yMean=mean(y);

SxSy= sum (x.*y) - n*xMean*yMean;
SySy= sum (y.*y) - n*yMean*yMean;
SxSx= sum (x.*x) - n*xMean*xMean;

CalcSlope=SxSy/SxSx;
CalcIntercept=yMean-CalcSlope*xMean;
StdDev=sqrt((SySy-CalcSlope*SxSy)/(n-2)); 
%StdDev=sqrt((SySy-(SxSy^2/SxSx))/(n-2))

fprintf(1,...
'LS Calc.\t%6.2f Slope\t%6.2f Intercept %6.2f Standard Deviation\n',...
    CalcSlope, CalcIntercept, StdDev);

kappa = -log(CalcSlope)/delta;
alpha = CalcIntercept/(1-CalcSlope);
sigma = StdDev * sqrt(-2*log(CalcSlope)/(delta*(1-CalcSlope^2)));

end

    
    
    
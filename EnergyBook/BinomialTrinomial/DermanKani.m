function DermanKani(ConstVol)
%DermanKani binomial fits price/strike/expiration 
%to tree based on Constant Volatility input parameter
%or fits local volatility to European Call Data
%example: XEO european option data of S&P100 index

close all
clc

set(0,'defaultaxeslinewidth',3); set(0,'defaultlinelinewidth',3);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold') 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

if (nargin == 1), %Check for Constant Volatility Input 
    ConstVolFlag=1  %when equal 1 use constant volatility  
else
    ConstVolFlag=0
end

r=0.02; %Libor

if (ConstVolFlag==1) %DK Code Can replicate CRR Tree
        Vol=ConstVol;
        d=0; %Assume No Dividends or Div. Yield
        S0=100;
        maxvol=1.1*ConstVol %Check later for stability
else  %Extract Implied Volatility from XEO data
    S0=567.10 %XEO Index Price as of 12/27/2010
    d=0.02; %Estimated dividend yield
    CallKvTdata=load('SP100European.dat'); %read directly 
    %Input Matrix should have 3 columns: time , strike, price
    Tvector=(CallKvTdata(:,1));
    Kvector=(CallKvTdata(:,2));
    Pvector=(CallKvTdata(:,3));
    CallFlagVector=(CallKvTdata(:,4));

    for i = 1: length(Tvector)
    VolVector(i)=InvertOptionPrice(r,S0, d, Tvector(i),...
        Kvector(i),Pvector(i),CallFlagVector(i));
    end
    maxvol=max(VolVector) %Check later for stability
    %%%%%begin extracted volatility plot
    %2D figure;%plot(Kvector,VolVector,'.','MarkerSize',15);%Tvector,
    figure  %3D
    Tlin=linspace(0,max(Tvector),51); %or min(Tvector)
    Klin=linspace(min(Kvector),max(Kvector),51);
    [Tgrid,Kgrid]=meshgrid(Tlin,Klin);
    %regularly spaced matric of data
    Volgrid=griddata(Tvector,Kvector,VolVector,Tgrid,Kgrid,'linear');

    mesh (Tgrid,Kgrid,Volgrid)
    axis tight; hold on

    plot3(Tvector,Kvector,VolVector,'.','MarkerSize',15);
    hold on

    Svector=S0*ones(1,length(VolVector));
    VolLine=min(VolVector)*ones(1,length(VolVector));
    %draw line at Current Price
    plot3(Tvector,Svector,VolLine,'green');

    xlabel('Expiration [Years]');
    ylabel('Strike'); 
    zlabel('Volatility Estimated');
    title('XEO S&P100 Out of the Money Option Data');
    %%%%%end extracted volatility plot
end  %%%if-else for volatility

%%%%%%%%%Implement Derman and Kani Algorithm%%%%%%%%%%%
Steps=5;
Tmax=2; %years
dt=Tmax/(Steps)
Prob=zeros(Steps,Steps); %pre-allocate
S=zeros(Steps+1,Steps+1);
Arrow=zeros(Steps+1,Steps+1);
S(1,1)=S0;
Arrow(1,1)=1;
erdt=exp(r*dt);
F(1,1)=S0*erdt;

for t = 1 : Steps;
    time=t*dt ;
    N=t; % equate steps in DK tree to steps used in binomial tree

%%%%%%Central Nodes%%%%%%%%%    
    if (rem(t+1,2)==0)    %current t level odd, 
  %next t+1 level = even # of nodes 
  %-> 2 node straddle center of tree
        midtop = (t+1)/2;
        K=S(midtop,t) ;
        if (ConstVolFlag==1)
            vol=ConstVol;
        else
            vol=max(0.1,...
                griddata(Tvector,Kvector,VolVector,time,K,'linear'))
        end
        
        CallFlag=1; EuropeanFlag=1; %European Call
        Call= VectorBinomial(S0,K,r,vol,time,N,d,CallFlag,EuropeanFlag);
   
        NodeSum=sum(Arrow(1:midtop-1,t)...
            .*(erdt.*S(1:midtop-1,t)-S(midtop,t)));
       
        F(midtop,t)=S(midtop,t)*erdt;
        S(midtop,t+1)=(S(midtop,t)...
           *(erdt*Call+Arrow(midtop,t)*S(midtop,t)-NodeSum))...
           / (Arrow(midtop,t)*F(midtop,t)-erdt*Call+NodeSum);

        midbot= midtop+1;
        S(midbot,t+1)=S0^2/S(midtop,t+1) ;

    else
        midtop = ceil((t+1)/2); %set bottom mode for call loop 
        midbot=midtop; %set upper limit for put loop
        S(midtop,t+1)=S0; %center trunk of tree        
    end

%%%Top of tree, Call option Nodes in the money%%%%%%%%
    for j=midtop-1:-1:1
        K=S(j,t); 
        if (ConstVolFlag==1)
            vol=ConstVol;
        else
            vol=max(0.1,...
                griddata(Tvector,Kvector,VolVector,time,K,'linear'));
        end
        
        CallFlag=1; EuropeanFlag=1; %European Call
        Call = VectorBinomial(S0,K,r,vol,time,N,d,...
                        CallFlag,EuropeanFlag);
        NodeSum=sum(Arrow(1:j-1,t).*(erdt.*S(1:j-1,t)-S(j,t)));
        F(j,t)=S(j,t)*erdt;
        S(j,t+1)=(S(j+1,t+1)*(erdt*Call-NodeSum)...
            -Arrow(j,t)*S(j,t)*(F(j,t)-S(j+1,t+1)))...
            /(erdt*Call-NodeSum-Arrow(j,t)*(F(j,t)-S(j+1,t+1)));
        if ( (S(j,t+1)<F(j,t)))...
                  % || S(j,t+1)> exp(maxvol*sqrt(dt))*F(j,t) ) 
            %Over-ride with Log spacing for stability
            %S(j,t+1)=S(j,t)^2/S(j+1,t+1);
            S(j,t+1)=exp(vol*sqrt(dt))*S(j,t);%More Stable
        end
    end%%% End Top Loop
    
%%%Bottom of tree, Put option Nodes in the money%%%%%%%%
    for j=(midbot+1):1:(t+1)
        K=S(j-1,t);
        if (ConstVolFlag==1)
            vol=ConstVol;
        else
            vol=max(0.1,...
                griddata(Tvector,Kvector,VolVector,time,K,'linear'));
        end
        CallFlag=0; EuropeanFlag=1; %European Put
        Put = VectorBinomial(S0,K,r,vol,time,N,d,...
                        CallFlag,EuropeanFlag);
        NodeSum=sum(Arrow(j:t,t).*(S(j-1,t)-erdt.*S(j:t,t)));
        F(j-1,t)=S(j-1,t)*erdt;
        S(j,t+1)=(S(j-1,t+1)*(erdt*Put-NodeSum)...
            +Arrow(j-1,t)*S(j-1,t)*(F(j-1,t)-S(j-1,t+1)))...
            /(erdt*Put-NodeSum+Arrow(j-1,t)*(F(j-1,t)-S(j-1,t+1)));
        if ( (S(j,t+1)>F(j-1,t)))%...
                %|| S(j,t+1)<exp(-maxvol*sqrt(dt))*F(j-1,t) )
            %Over-ride with Log spacing to avoid neg prob
               % S(j,t+1)=S(j-1,t)^2/S(j-1,t+1) ;
               S(j,t+1)=exp(-vol*sqrt(dt))*S(j-1,t);%More Stable
        end
    end%%% End Bottom Loop
  
 %%%%Arrow Debreu Tree and Node Probability
  %Vectorize -> Calculate all probabilities at this level   
     Prob(1:t,t)=(F(1:t,t)-S(2:t+1,t+1))./(S(1:t,t+1)-S(2:t+1,t+1));
  %Calculate Top and Bottom Arrow Debreu nodes for next time step
     Arrow(1,t+1)=Arrow(1,t)*Prob(1,t)/erdt; %up
     Arrow(t+1,t+1)=Arrow(t,t)*(1-Prob(t,t))/erdt; %down
  %Vectorize -> Calculate rest of Arrow Debreu at next time step 
 	 Arrow(2:t,t+1)=(Arrow(1:t-1,t).*(1-Prob(1:t-1,t))...
         +Arrow(2:t,t).*Prob(2:t,t))/erdt;
end %%%% End Derman Kani Calculation Time loop
  
figure %
for t=1:Steps
    for i=1:t
       clear x y
       x(:,1)=[(t-1)*dt t*dt];%+dt
       x(:,2)=x(:,1); 
       y(:,1)=[S(i,t) S(i, t+1)];
       y(:,2)=[S(i,t) S(i+1,t+1)];
       plot(x,y)
       txstr1(1) = {num2str(S(i,t),'S=%.1f')};
       txstr1(2) = {num2str(Arrow(i,t),'AD=%.2f')}; 
       txstr1(3) = {num2str(Prob(i,t),'p=%.2f')}; 
       %txstr1(4) = {num2str(F(i,t),'%.2f')}; 
       text((t-1)*dt,S(i,t),txstr1);
       hold on
    end
end
t=Steps+1;
    for i=1:t
     % Expiration nodes
       txstr2(1) = {num2str(S(i,t),'S=%.1f')};
       txstr2(2) = {num2str(Arrow(i,t),'AD=%.2f')}; 
       text((t-1)*dt*0.95,S(i,t),txstr2);
       hold on
    end
axis tight  
xlabel('Time (Years)'); ylabel('Asset Value');
hold off
end  % DermanKani function

function Vol = InvertOptionPrice(r,S,d,T,K,P,CallFlag);
    VolGuess=0.2;
    options=[];%=optimset('Display','iter');
    if (CallFlag==1) %Call
        [Estimates]=fminsearch(@callError,VolGuess,...
            options,K,S,T,r,d,P);
    else %(CallFlag==0) %Put
        [Estimates]=fminsearch(@putError,VolGuess,...
            options,K,S,T,r,d,P);
    end
    Vol=max(0.001,Estimates);
end
   
function error = callError (vol,K,S0,T,r,d,P);
    N=50; CallFlag=1; EuropeanFlag=1;
    c= VectorBinomial(S0,K,r,vol,T,N,d,...
                        CallFlag,EuropeanFlag);    
    error=(c-P)^2;% minimize the squared error
end

function error = putError (vol,K,S0,T,r,d,P);
    N=50; CallFlag=0; EuropeanFlag=1;
    p= VectorBinomial(S0,K,r,vol,T,N,d,...
                        CallFlag,EuropeanFlag);
    error=(p-P)^2;% minimize the squared error
end



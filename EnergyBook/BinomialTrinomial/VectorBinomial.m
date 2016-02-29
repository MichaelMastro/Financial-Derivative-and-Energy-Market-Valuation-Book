function Val=VectorBinomial(S0,K,r,sigma,T,N,d,...
                        CallFlag,EuropeanFlag)
%VectorBinomial computes option price via binomial tree
%for Call (CallFlag=1) or Put (CallFlag=0) 
%for European (EuropeanFlag=1) or American (EuropeanFlag=0) 
%Vector of asset (S) and option (V) prices greatly reduces
%CPU and Memory Constraints compared to Array for large
%number of steps (N)
%A N step binomial tree has 2N+1 nodes with only
%N distinct asset prices possible

if (nargin == 0), %Check for Data Input 
  S0=50;
  K=50;
  r=0.1;
  sigma=0.5;
  T=1;
  N=100;
  d=0;
  CallFlag=1;
  EuropeanFlag=0;
end

dt=T/N;
a=exp((r-d)*dt);
u=exp(sigma*sqrt(dt)); %standard dev = annual vol. x sqrt(dt)
d=exp(-sigma*sqrt(dt));
p=(a-d)/(u-d);
S=zeros(2*N+1,1);
V=zeros(2*N+1,1);
DiscUpProb=exp(-r*dt)*p;
DiscDownProb=exp(-r*dt)*(1-p);
%Initial Asset price and center trunk of tree
S(N+1)=S0;

downpow=[N:-1:1];;%Vectorize down movements
uppow=[1:N];%Vectorize up movements
%S(1) is lowest price S0*d*d*d*d...
S(1:N)=d.^downpow*S0;
%S(2N+1) is highest price S0*u*u*u...
S(N+2:2*N+1)=u.^uppow*S0;
%Node 1,3,...2N+1 are expiration nodes
%Value at expiration is excess over Strike Price K
if (CallFlag==1)
    V(1:2:2*N+1)=max(0,S(1:2:2*N+1)-K);
else
    V(1:2:2*N+1)=max(0,K-S(1:2:2*N+1)); 
end

if (EuropeanFlag==1) %Propagate European Call or Put 
    for Step=1:N
    V(Step+1:2:(2*N+1-Step))...
        =DiscUpProb*V(Step+2:2:(2*N+2-Step))...
            +DiscDownProb*V(Step:2:(2*N-Step));
    end
elseif (CallFlag==1) %American Call
    for Step=1:N
    V(Step+1:2:(2*N+1-Step))...
        =max(S(Step+1:2:(2*N+1-Step))-K,...
            (DiscUpProb*V(Step+2:2:(2*N+2-Step))...
            +DiscDownProb*V(Step:2:(2*N-Step))));
    end
else        %American Put
for Step=1:N
    V(Step+1:2:(2*N+1-Step))...
        =max(K-S(Step+1:2:(2*N+1-Step)),...
            (DiscUpProb*V(Step+2:2:(2*N+2-Step))...
            +DiscDownProb*V(Step:2:(2*N-Step))));
end
end
Val=V(N+1);
end



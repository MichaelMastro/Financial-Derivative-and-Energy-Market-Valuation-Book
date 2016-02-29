function V=VectorTrinomial(S0,K,r,sigma,T,N,d,...
                        CallFlag,EuropeanFlag)
%VectorTrinomial computes option price via trinomial tree
%for Call (CallFlag=1) or Put (CallFlag=0) 
%for European (EuropeanFlag=1) or American (EuropeanFlag=0) 
%Vector of asset (S) and option (V) prices greatly reduces
%CPU and Memory Constraints compared to Array for large
%number of steps (N)
%Each new time step in the tree is identical to the 
%previous time except one higher and one lower node is added
%N step trinomial tree has 2N+1 nodes at expiration

clc
close all

if (nargin == 0), %Check for Data Input 
  S0=50;
  K=50;
  r=0.1;
  sigma=0.5;
  T=1;

  d=0; %dividend yield
  CallFlag=0;
  EuropeanFlag=1; 
end

X0=log(S0);

dt=T/N;
dx=sigma*sqrt(3*dt);
nu=r-0.5*sigma^2-d;

X=zeros(2*N+1,1); %Asset
V=zeros(2*N+1,1); %Option  

j=-N:1:N;
X(j+N+1)=X0-j*dx;
S=exp(X);
if (CallFlag==1)
    V=max(S-K,0);
else
    V=max(K-S,0);
end

%Only use probabilities when discounting option value to 
%give present value at previous time step
DiscPu=exp(-r*dt)*0.5*( (nu^2*dt^2+sigma^2*dt)/dx^2 + nu*dt/dx);
DiscPd=exp(-r*dt)*0.5*( (nu^2*dt^2+sigma^2*dt)/dx^2 - nu*dt/dx);
DiscPm=exp(-r*dt)*(1-( (nu^2*dt^2+sigma^2*dt)/dx^2));

if (EuropeanFlag==1)
    for edge=0:N-1
        up=(1:1: 2*N+1 -2*edge -2); %Vup=pu*V(up);
        %mid=(up+1); Vmid=pm*V(mid);
        %down=(up+2); Vdown=pd*V(down);
        %Vprev=Vup+Vmid+Vdown
        Vprev=DiscPu*V(up)+DiscPm*V(up+1)+DiscPd*V(up+2);
        V=Vprev;
    end
elseif (CallFlag==1) %American Call         
    for edge=0:N-1
        up=(1:1: 2*N+1 -2*edge -2); 
        current=(edge+2 :1: 2*N+1-1-edge);
        Vprev=DiscPu*V(up)+DiscPm*V(up+1)+DiscPd*V(up+2);
        V=max(S(current)-K,Vprev);
    end
else %American Put
    for edge=0:N-1
        up=(1:1: 2*N+1 -2*edge -2); 
        current=(edge+2 :1: 2*N+1-1-edge);
        Vprev=DiscPu*V(up)+DiscPm*V(up+1)+DiscPd*V(up+2);
        V=max(K-S(current),Vprev);
    end    
end
end




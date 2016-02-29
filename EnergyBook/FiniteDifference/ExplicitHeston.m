function [c0  C0] = ExplicitHeston(S0,K0,r,dt,T,...
    kappa,lambda,sigma,rho,eta,v0)
% ExplicitHeston: explicit finite difference
% of Heston Model. Two-dimensional grid of coupled
% stochastic volatility and asset price
% Time step backwards from known call option value 
% at maturity to time zero. Explicit FD requires fine
% time step grid. Easier to let code calculate dt.

close all;
clc;
set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold'); 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

disp ('Heston Stochastic Volatility Option Calculation')
disp ('via Explicit Finite Difference')

if (nargin < 10), v0 = 0.5; end
if (nargin < 9), eta = 0.4; end
if (nargin < 8), rho = -0.8; end
if (nargin < 7), sigma = 1.5; end
if (nargin < 6), lambda =0.05; end 
if (nargin < 5), kappa = 1; end
if (nargin < 4), T = 1, end
if (nargin < 4), dt = T/10; end %Recalculate Below
if (nargin < 3), r = 0.03; end
if (nargin < 2), K0 = 50, end
if (nargin < 1), S0 = 50; end

str1=['\kappa=' num2str(kappa)];
str2=[' \sigma=' num2str(sigma)];
str3=[' \lambda=' num2str(lambda)];
str4=[' \eta=' num2str(eta)];
str5=[' \rho=' num2str(rho)];
textstr=[str1 str2 str3 str4 str5];

I=30 %I+1 = number of Stock Index Points
J=30 %J+1 = number of Volatility Index Points

Smax=2*K0
dS=Smax/I
S=[0:I]*dS;

vMax=1.5*(v0+sigma+eta);
dv=vMax/J
v=[0:J]*dv;

% Check if dt less than stability level
% if not stable then change dt
stab=1/(I^2*vMax + J*sigma^2+r)
if (dt<=stab)
    fprintf('Stable: [dt = %f] < ',dt)
    fprintf('[1/(I^2*vMax + J*sigma^2+r) = %f]\n',stab)
else
    fprintf('Not Stable: [dt = %f] > ',dt)
    fprintf('[1/(I^2*vMax + J*sigma^2+r) = %f]\n',stab)
    dt=stab/2;
    fprintf('Set dt = (stability requirement)/2 = %f\n',dt)
end

Tsteps=round(T/dt)+1; %e.g., T=1, dt=1 -> Tsteps=2
jInt=2:J;   %Interior Volatility Index Points
iInt=2:I;   %Interior Stock Index Points
c1=zeros(I+1,J+1); c2=zeros(I+1,J+1); % European
%C1=zeros(I+1,J+1); C2=zeros(I+1,J+1); % American
SK=repmat(max(S-K0,0)',1,J+1);
c2=SK; %C2=c2;
for i=iInt
 for j=jInt
% Calculate Coefficients used in Time Stepping  
  A(i,j)=-i^2*v(j)*dt -sigma^2*j*dt/dv +1 -r*dt;
  C(i,j)=(0.5*i^2*v(j) - 0.5*r*i)*dt;
  D(i,j)=(0.5*i^2*v(j) + 0.5*r*i)*dt;
  E(i,j)=(0.5*sigma^2*j/dv -...
      0.5*(kappa*(eta-v(j))-lambda)/dv)*dt;
  F(i,j)=(0.5*sigma^2*j/dv +....
      0.5*(kappa*(eta-v(j))-lambda)/dv)*dt;
  B(i,j)=0.25*rho*sigma*i*j*dt;
 end
end

disp('Time Step')
sqrtTsteps=sqrt(Tsteps);
% Time Steps Backwards from known Intrinsic Value at
% Expiration
tic
for n=Tsteps:-1:1
    if (mod(n,sqrtTsteps)<0.5)
        disp(n) %crude partial time counter display 
    end 
 for i=iInt
  for j=jInt
c1(i,j)=A(i,j)*c2(i,j) +C(i,j)*c2(i-1,j) +D(i,j)*c2(i+1,j)...
  +E(i,j)*c2(i,j-1) +F(i,j)*c2(i,j+1)...
  +B(i,j).*(c2(i-1,j-1)-c2(i-1,j+1)-c2(i+1,j-1)+c2(i+1,j+1));
%C1(i,j)=A(i,j)*C2(i,j) +C(i,j)*C2(i-1,j) +D(i,j)*C2(i+1,j)...
 % +E(i,j)*C2(i,j-1) +F(i,j)*C2(i,j+1)...
 % +B(i,j).*(C2(i-1,j-1)-C2(i-1,j+1)-C2(i+1,j-1)+C2(i+1,j+1));
% Check Early Exercise of American Call (C1)
%C1(i,j)=max(SK(i,j),C1(i,j)); 
  end
 end
% BC at S=0
    c1(1,1:J+1)=0; %C1(1,1:J+1)=0; 
% BC at S=Smax dc/ds=1
    c1(I+1,1:J+1)=c2(I,1:J+1)+dS; 
    %C1(I+1,1:J+1)=C2(I,1:J+1)+dS; 
% BC at v=vMax  dc/dv=0
    c1(1:I+1,J+1)=c2(1:I+1,J);
   % C1(1:I+1,J+1)=C2(1:I+1,J);
% BC at v=0  (dc/dt)+kappa*nu*(dc/dv)=0
c1(1:I+1,1)=(1-kappa*eta*dt/dv)*c2(1:I+1,1)...
              + kappa*eta*dt/dv*c2(1:I+1,2) ;
%C1(1:I+1,1)=(1-kappa*eta*dt/dv)*C2(1:I+1,1)...
 %             + kappa*eta*dt/dv*C2(1:I+1,2) ;
c2=c1; %C2=C1;  % Set for Next Time Step
end
explicitTime=toc;
fprintf('Explicit Heston Calculation Time: %f\n',...
    explicitTime);
% Interpolation of European (c0) and American (C0) call
% values at S0 and v0
c0 = interp2 (S, v, c1, S0, v0);
%C0 = interp2 (S, v, c1, S0, v0);

PlotFlag=1; % Optionally plot results
if (PlotFlag==1)
    Smatrix=repmat(S',1,J+1);
    vMatrix=repmat(v,I+1,1);
    figure; 
    surf (Smatrix,vMatrix,c2)
    STstr = ['Strike = ' num2str(K0) ', T = ' num2str(T)];
    title(['Heston Finitie Difference: ' STstr])
    zlabel ('European Call'); xlabel ('Asset Price'); 
    ylabel ('Variance'); axis tight
    text(0,vMatrix(I,J)*1.3,c2(I,J)*1.3,textstr)
end
end
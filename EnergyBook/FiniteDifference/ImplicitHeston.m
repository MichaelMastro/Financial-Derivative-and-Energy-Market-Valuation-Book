function [c0 ImplicitTime] = ImplicitHeston(S0,K0,r,dt,T,...
    kappa,lambda,sigma,rho,eta,v0)
% ImplicitHeston: implicit finite difference
% of Heston Model. Two-dimensional grid of coupled
% stochastic volatility and asset price
% Time step backwards from known call option value 
% at maturity to time zero. Implicit FD relaxes 
% time step requirement that slows explicit scheme
% Implicit scheme calculates several unknown points
% in current time layer consistently while using only
% one known (forward in time) point
% Handling the boundary conditions is difficult 
% in a matrix based calculation of an implicit or 
% Crank-Nicolson representation of two-dimensional 
% Heston or similar 2D model
% To greatly simplify the implementation this code
% uses a successive over relation (SOR) approach. 

%close all;
%clc;
set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold'); 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

disp ('Heston Stochastic Volatility Option Calculation')
disp ('via Implicit Finite Difference')

if (nargin < 11), v0 = 0.5; end
if (nargin < 10), eta = 0.4; end
if (nargin < 9), rho = -0.8; end
if (nargin < 8), sigma = 1.5; end
if (nargin < 7), lambda =0.05; end 
if (nargin < 6), kappa = 1; end
if (nargin < 5), T = 1, end

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

if (nargin < 4) % Calculate dt
    stab=1/(I^2*vMax + J*sigma^2+r)
    fprintf('[1/(I^2*vMax + J*sigma^2+r) = %f]\n',stab)
    % Implicit allows large multiple of stability level
    dt=stab*20
end 

Tsteps=round(T/dt)+1; %e.g., T=1, dt=1 -> Tsteps=2
jInt=2:J;   %Interior Volatility Index Points
iInt=2:I;   %Interior Stock Index Points
c1new=zeros(I+1,J+1); c1old=zeros(I+1,J+1);
c2=zeros(I+1,J+1); % European
SK=repmat(max(S-K0,0)',1,J+1);
c2=SK; %C2=c2;
for i=iInt
 for j=jInt
% Calculate Coefficients used in Time Stepping
% Implicit coeffecients are equal to
% negative of coefficients in explicit scheme
% except a=Ap2=-A+2
  nAp2(i,j)=1+i^2*v(j)*dt +sigma^2*j*dt/dv +r*dt;
  nC(i,j)=-(0.5*i^2*v(j) - 0.5*r*i)*dt;
  nD(i,j)=-(0.5*i^2*v(j) + 0.5*r*i)*dt;
  nE(i,j)=-(0.5*sigma^2*j/dv -...
      0.5*(kappa*(eta-v(j))-lambda)/dv)*dt;
  nF(i,j)=-(0.5*sigma^2*j/dv +....
      0.5*(kappa*(eta-v(j))-lambda)/dv)*dt;
  nB(i,j)=-0.25*rho*sigma*i*j*dt;
 end
end

c1new=c2; %Initialize with time T terminal values
counter = 0; err=1000;
disp('Time Step   Counter  Error')
sqrtTsteps=sqrt(Tsteps);
% Time Steps Backwards from known Intrinsic Value at
tic
% Expiration
for n=Tsteps:-1:1
    if (mod(n,sqrtTsteps)<0.5)
        ([n counter err]) %crude partial time counter display 
    end

% BC at S=0
    c1new(1,1:J+1)=0; %C1(1,1:J+1)=0; 
% BC at S=Smax dc/ds=1
    c1new(I+1,1:J+1)=c2(I,1:J+1)+dS; 
    %C1(I+1,1:J+1)=C2(I,1:J+1)+dS; 
% BC at v=vMax  dc/dv=0
    c1new(1:I+1,J+1)=c2(1:I+1,J);
   % C1(1:I+1,J+1)=C2(1:I+1,J);
% BC at v=0  (dc/dt)+kappa*nu*(dc/dv)=0
    c1new(1:I+1,1)=(1-kappa*eta*dt/dv)*c2(1:I+1,1)...
              + kappa*eta*dt/dv*c2(1:I+1,2) ;
%C1(1:I+1,1)=(1-kappa*eta*dt/dv)*C2(1:I+1,1)...
 %  + kappa*eta*dt/dv*C2(1:I+1,2) ;  
 
 counter = 0; err=10000;
 c1old=c2; 
 w=1.3; %Relaxation parameter
 % SOR or Gauss-Seidel (w=0) iteratively approaches
 % true values for c in M*c1=c2 where c2 is the 
 % previously calculation (forward in time) layer
 % Just calculated (updated) values c1new are used
 % else previous iteration value c1old are used
 while ( (counter < 100) && (err >0.01))
     for i=iInt
      for j=jInt
         c1new(i,j)= (1-w)*c1new(i,j)+...
                (w/nAp2(i,j))*(c2(i,j)...
            -nC(i,j)*c1new(i-1,j) -nD(i,j)*c1old(i+1,j)...
            -nE(i,j)*c1new(i,j-1) -nF(i,j)*c1old(i,j+1)...
            -nB(i,j).*(c1new(i-1,j-1)-c1new(i-1,j+1)...
            -c1old(i+1,j-1)+c1old(i+1,j+1)));
      end
     end
    err=norm(c1old-c1new);
    counter=counter+1;
    c1old=c1new;
    %disp ([n counter err])
 end
 c2=c1new;% Set for Next Time Step
end
ImplicitTime=toc;
fprintf('Implicit Calculation Time: %f\n',...
    ImplicitTime); 
% Interpolation of European (c0)  call
% values at S0 and v0
c0 = interp2 (S, v, c1new, S0, v0);

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
function [C,k]=FFToptionOTM (CharFunc, Param, alpha, wEnd, N)
% FFToption Performs Carr and Madam Call FFT Option 
% Calculation. FFToption returns array of Call Prices
% and their corresponding strikes for a fixed S0.
% Calculation 

if (nargin < 3), alpha=1.25; end
if (nargin <4),wEnd=500; end
if (nargin <5),N=2^12; end

% characteristic function discretization 
% Increase Max Omega when more kurtosis
% frequency grid size for char func
w=linspace(0,wEnd,N);
dw=w(2)-w(1)   
dk=2*pi/(dw*N)
k1=-(dk*(N)/2)  
kend=(dk*(N-1)/2)   
k=k1:dk:kend; %linspace(k1,kend,N)
klength=length(k)

r=Param.r;
T=Param.T;

wN=w-i*alpha;
wP=w+i*alpha;

PsiN=exp(-r*T).*( (1./(1+i*wN)) -...
    (feval(CharFunc, (-i), Param)./(i*wN)) -...
    ( feval(CharFunc, (wN-i), Param) ./ (wN.^2 -i*wN) ) );

PsiP=exp(-r*T).*( (1./(1+i*wP)) -...
    (feval(CharFunc, (-i), Param)./(i*wP)) -...
    ( feval(CharFunc, (wP-i), Param) ./ (wP.^2 -i*wP) ) );

Psi = (PsiN - PsiP) /2;

%Trapezoial
% q=ones(1,N); q(1)=0.5; q(end)=0.5; 
% Carr and Madan suggest Simpson's Quadratic weight
% as Euler terms can be highly Oscillatory
q=3+(-1).^(1:N); q(1)=1; q=q/3;

FFTinput=exp(-i*w*k1).*q.*Psi*dw;

%C=(exp(-alpha*k) / pi) .* real(fft(FFTinput));

C=(1./(pi*sinh(alpha*k))).* real(fft(FFTinput));
end
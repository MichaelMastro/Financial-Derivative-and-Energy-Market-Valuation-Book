function [C,k]=FFToption(CharFunc, Param, alpha, wEnd , N)
% FFToption Performs Carr and Madam Call FFT Option 
% Calculation. FFToption returns array of Call Prices
% and their corresponding strikes for a fixed S0.
% Calculation 

if (nargin < 3), alpha=1.25; end

if (nargin <4),wEnd=500; end

if (nargin <5), N=2^13; end

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

Psi=feval(CharFunc, (w-(alpha+1)*i), Param);
% 'Spiral' Oscillation in complex Fourier space 
%figure
%subplot (2,1,1); 
%plot3 (w, real(Psi), imag(Psi))

Psi=Psi.*exp(-r*T)./...
    ( alpha^2+alpha-w.^2 + i*(2*alpha+1).*w);

%Trapezoial
% q=ones(1,N); q(1)=0.5; q(end)=0.5; 
% Carr and Madan suggest Simpson's Quadratic weight
% as Euler terms can be high
q=3+(-1).^(1:N); q(1)=1; q=q/3;

FFTinput=exp(-i*w*k1).*q.*Psi*dw;
% FFTinput much more complex than output of 
% Characteristic function. Larger N reveals much 
% more detail
%subplot (2,1,2);
%plot3 (w, real(FFTinput), imag(FFTinput))
%xlim ([0, 10])

if (nargin <5)
    C=(exp(-alpha*k) / pi) .* real(fft(FFTinput));
else
    C=(exp(-alpha*k) / pi) .* real(frfft(FFTinput,1/N));
end
end




function [CDF FFTinput omega x] =...
    cdfCharFuncFFT (p, N, OmegaEnd,eta); 
% CharFuncInv produces CDF via FFT of
% alpha-stable characteristic function
% Decay factor eta necessary to avoid divergence at zero

if (nargin <1)
p.t = 1/12;
p.mu=0;
p.alpha=2;
p.beta=0;
p.sigma =0.3;
end

if (nargin <2), N=2^9
end
if (nargin <3),OmegaEnd=100
end
if (nargin <4), eta=1
end

% characteristic function discretization
% Increase Max Omega when more kurtosis
% frequency grid size for char func
omega=linspace(0,OmegaEnd,N)';
DeltaOmega=omega(2)-omega(1) ;   

dx=2*pi/(DeltaOmega*N);
xpoint=(dx*(N-1)/2);
x=(-xpoint:dx:xpoint);

ShiftOmega=eta+i*omega+eps*i;

CDFphi=((1+eps*i)./(eta-i*omega)).*phi_alpha(ShiftOmega,p);
FFTinput=exp(-i*x(1)*omega).*CDFphi*DeltaOmega;
% Trapezoidal Rule: 1st and last divided by 2
FFTinput(1)=FFTinput(1)/2; FFTinput(end)=FFTinput(end)/2;
CDF=(exp(eta*x)/pi)'.*real(fft(FFTinput));
CDF=(CDF-min(CDF))./max(CDF); %Normalize
end


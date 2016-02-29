function [PDF FFTinput omega x] = CharFuncFFT (p, N, OmegaEnd); 
%CharFuncInv produces PDF via FFT of
%alpha-stable characteristic function

if (nargin <1)
p.t = 1;
p.mu=0.1;
p.alpha=1.9;
p.beta=0.1;
p.sigma =0.1;
end

if (nargin <2), N=2^9
end
if (nargin <3),OmegaEnd=100
end


%characteristic function discretization
%Increase Max Omega when more kurtosis
%frequency grid size for char func
omega=linspace(0,OmegaEnd,N)';
DeltaOmega=omega(2)-omega(1) ;   

dx=2*pi/(DeltaOmega*N);
xpoint=(dx*(N-1)/2);
x=(-xpoint:dx:xpoint);

phia=phi_alpha(omega,p);
FFTinput=exp(-i*x(1)*omega).*phia*DeltaOmega;
%Trapezoidal Rule: 1st and last divided by 2
FFTinput(1)=FFTinput(1)/2; FFTinput(end)=FFTinput(end)/2;
PDF=(1/pi)*real(fft(FFTinput));

end






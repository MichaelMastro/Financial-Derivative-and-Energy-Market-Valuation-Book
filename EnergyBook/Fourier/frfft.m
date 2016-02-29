function y = frfft (X, alpha)
% frfft fractional fast Fourier transform called via
% Chirp z transform function mczt
% Defaults to standard FFT when alpha=1/N 

    if (nargin <2), alpha=1./N; end

    N=length(X); 
    %set number of output bins M = number of input bins N

    W=exp(-i*2*pi*alpha);%exp(-i*2*pi*phi);
    % i.e., decrease alpha -> decrease frequency range
    
    A=1; %exp(i*2*pi*theta)
    method ='Bluestein'; 
    y = mczt ( X, N , W, A, method); % Call Chirp z (M=N)
end

    
function [y] = mczt ( X, M , W, A, method);
% mczt calculates chirp z transform similar to czt in
% Matlab's Signal Processing Toolbox
% 'Direct' and 'convolution' methods are for
% education. 'If' statements slow down code 

if (nargin < 1)  % Use for Test
    N=2^6; n=(0:N-1);
    f1=sqrt(N); f2=1.1*f1;
    X=sin(f1*2*pi*n/N)+sin(f2*2*pi*n/N);
else
    N=length (X); n=(0:N-1);
end

if (nargin <2), M=N; end % Output length
M=round(M);
Lexp = nextpow2(N+M); L=2^Lexp;

if (nargin <3), W=exp(-i*2*pi/M); 
end %Default FFT if M=N

if (nargin <4), A=1; end %=exp(0); end
if (nargin <5), method = 'Bluestein';   end

% direct method should be faster for very small
% output arrays but Matlab FFT routines are much
% fast than this crude code
if (strcmpi(method , 'direct'))
    Wvector = W.^n;
    Avector= A.^n;
    for k=1:M
        km1=k-1;
        y(k)=sum(X.*Avector.*Wvector.^(-km1));
    end
    
elseif (strcmpi(method , 'convolution')) %Slowest
    WvectorSQ = W.^((n.^2)/2);
    Avector= A.^(-n);
    for k=1:M
        km1=k-1;
        y(k)=W.^((km1.^2)/2)...
            .*sum(X.*Avector.*WvectorSQ...
                .*W.^((-(km1-n).^2)/2));
    end
    
else % 'Bluestein' = Bluestein iFFT / 2-FTT process
    WvectorSQ = W.^((n.^2)/2);
    Avector= A.^(-n);
    x(1,:)=X(:); % make dimension same as A and W

    z1=[ x.*Avector.*WvectorSQ   zeros(1,L-N) ];
    
    LNto1=(L-N):-1:1;
    WvectorSQbackwards=W.^(-(LNto1.^2)/2);
    % First component in z2 is increased to length N
    % to fill 'arbitrary region 
    z2=[ WvectorSQ.^(-1)     WvectorSQbackwards ];
    
    fz1=fft(z1);
    fz2=fft(z2);
    fz=fz1.*fz2;
    ifz=ifft(fz);
    
    k=[0:(M-1)];
    WvectorSQk = W.^((k.^2)/2);
    y=WvectorSQk.*ifz(1:M);
end

end
function Y = GenAlphaDist(alpha,beta,mu,sigma,N);
%GenAlphaDist generates an array of length N drawn from
%an alpha stable distribution defined by 
%center: mu, scale: sigma; skew: beta, tail param: alpha

%Sz. Borak, W. Härdle, R. Weron (2005)'Stable distributions' 
%in 'Statistical Tools for Finance and Insurance', 
%eds. P. Cizek, W. Härdle, R. Weron, Springer-Verlag, 
%Berlin, 21-44; also see more robust 
%code by Weron http://fmwww.bc.edu/repec/bocode/s/stablernd.m

if (nargin<5), N=2^12; end
if (nargin<4), sigma=1; end
if (nargin<3), mu=0;   end
if (nargin<2), beta=0; end
if (nargin<1), alpha=2;end

%W = independent exponential random variable, 
%pdf=lambda*exp(-lambda*W) , mean=(1/lambda) 
%if mean =1 -> W =-log(n), n~Uniform(0,1)
W=-log(rand(1,N));
%random variable uniformly distributed on (-pi/2,pi/2)
U=pi*(rand(1,N)-0.5); 

if ((alpha>0.99) & (alpha<1.01)) %(alpha == 1)
    % alpha equal to one
    xi=pi/2; %lower case greek letter 'xi'
    X=(1/xi).*( (pi/2+beta.*U).*tan(U)...
        -beta.*log( ((pi/2).*W.*cos(U))./(pi/2+beta.*U)));
    Y=sigma*X+(2/pi).*beta.*sigma.*log(sigma)+mu;
else % alpha not equal to one
    zeta=-beta*tan(pi*alpha/2);
    xi=(1/alpha).*atan(-zeta);
    X=(1+zeta.^2).^(1/(2*alpha))...
        .*(sin(alpha.*(U+xi)) ./ (cos(U)).^(1/alpha))...
        .*(cos(U-alpha.*(U+xi))./W).^((1-alpha)/alpha);
    Y=sigma*X+mu; 
end
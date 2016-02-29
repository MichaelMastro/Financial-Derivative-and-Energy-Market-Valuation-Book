function [ C ] = OptionIntegration ( S0, K, r, T );
% OptionIntegration computes Call Option by numerical 
% integration of the characteristic equation of 
% Log-Normal (Gaussian log-return) model and Merton 
% Jump-Diffusion model. Results are compared to 
% analytical Black-Scholes (Gaussian log-return)  
% equation and Merton Series Expansion result.
% Numerical Integration technique is known to be
% 'stable' only for limited range of parameters
    clc
% Pass in or default to the 'important' 
    if (nargin <4), T=1; end
    if (nargin <3), r=0.05; end
    if (nargin <2), K=30; end 
    if (nargin <1), S0=50; end 
% For compactness, model input parameters are put into 
% InPar Structure
    InPar.volBS=0.9;   volBS=InPar.volBS;
    InPar.d=0;         d=InPar.d;

%%% Black-Scholes Analytical %%%
disp ('Black Scholes Model')
% BS = Current Stock price x probability [PI1(d1)=CDF(d1)] 
%      + discounted exercise payment x [PI2(d2)=CDF(d2)] 
    [cBS PI1BS PI2BS] = BlackScholesCall (K,S0,T,volBS,r,d);
% In Black-Scholes equation 
% PI1 = Delta (dC/dS) Hedge ratio for a Call
%   i.e, the rate of change in Call value for a change in S
% PI2 = Risk-Adjusted Probability that Option Exercised
%   i.e., finish in the money
    sprintf('Black-Scholes: PI1 %g , PI2 = %g', PI1BS,PI2BS) 
    sprintf('Black-Scholes: Call Price = %g', cBS ) 
%%% Black-Scholes Numerical %%%
% Use common PI function twice for PI1 (Number = 1)
% and PI2 (Number = 2)
    Model='Gaussian'; % Gaussian Log-Return Model
    CBSnumerical=S0*PI(S0,K,r,T, Model, 1, InPar)...
        -K*exp(-r*T)*PI(S0,K,r,T, Model, 2, InPar);
sprintf('Bakshi and Madam: Fourier Integration %s model = %g',...
                Model, CBSnumerical )

   CBSnumericalBates= S0 - exp(-r*T)*K*...
                (0.5+(1/pi)*Bates(S0,K,r,T, Model, 1, InPar));
   sprintf('Bates: Fourier Integration with %s model = %g',...
                Model, CBSnumericalBates )
            
%%% Merton Series %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp ('Merton Model')
    muJ=0.1            
    sigmaJ=0.1           
    lambda=0.1     
    cMertonSeries = MertonSeriesCall (K,S0,T,volBS,...
                        r, d, muJ, sigmaJ, lambda);
    sprintf('Merton Series: Call Price = %g', cMertonSeries) 

%%% Merton Numerical %%%
% Merton model input parameters put into InPar Structure
    InPar.muJ=muJ;
    InPar.sigmaJ=sigmaJ;
    InPar.lambda=lambda;
    Model='Merton';
    CMertonNumerical=S0*PI(S0,K,r,T, Model, 1, InPar)...
        -K*exp(-r*T)*PI(S0,K,r,T, Model, 2, InPar);
    sprintf('Bakshi and Madam Fourier Integration with %s model = %g',...
    Model, CMertonNumerical )

   CBSnumericalBates= S0 - exp(-r*T)*K*...
                (0.5+(1/pi)*Bates(S0,K,r,T, Model, 1, InPar));
   sprintf('Bates: Fourier Integration with %s model = %g',...
                Model, CBSnumericalBates )
end

function BatesVal = Bates ( S0, K, r, T, Model, Number, InPar)
% Bates performs integration via Matlab's quad or quadl function
% one integration and faster convergence with quadratic denominator
    BatesVal= quadl(@IntBates, realmin , 80,...
                [],[],S0,K,r,T, Model, Number, InPar );
% Bottom of intergration range is eps to avoid div by 0
end       

function IntBatesVal = IntBates (w,S0,K,r,T,Model,Number,InPar)
% IntFunction is the function 'Integrated' by the PI function 
% First, set phi to name of characteristic funciton 
    if (strcmp(Model,'Gaussian'))
        phi='phiBS';
        param.sigmaBS=InPar.volBS; 
        param.mu=r-0.5*param.sigmaBS^2; %risk-neutral
        param.T=T;
    elseif (strcmp(Model,'Merton'))
        phi='phiMerton';
        param.sigmaBS=InPar.volBS; sigmaBS=param.sigmaBS;
        param.muJ=InPar.muJ;       muJ=param.muJ;
        param.sigmaJ=InPar.sigmaJ; sigmaJ=param.sigmaJ;
        param.lambda=InPar.lambda; lambda=param.lambda;
% risk-neutral drift
        param.mu=r-0.5*sigmaBS^2 ...
                    -lambda*(exp(muJ+0.5*sigmaJ^2)-1); 
        param.T=T;
    else
        disp('Improper Model')
    end
    k=log(K/S0); % dimensionless moneyness
% Similar to Black-Scholes, PI1 calculated with stock as
% numeraire and PI2 calcuted with riskless bond as numeraire
        IntBatesVal=real(exp(-i*w*k).*feval(phi,w,param)...
                               ./(i.*w.*(1-i*w)));
end

function [ PIval ] = PI ( S0, K, r, T, Model, Number, InPar)
% PI performs integration via Matlab's quad or quadl function
% PI1 and PI2 are similar to Black-Scholes PI1=CDF(d1) and
% PI2 = CDF(d2)
    PIval= 0.5 + (1/pi)*quadl(@IntFunction, realmin, 80,...
                [],[],S0,K,r,T, Model, Number, InPar );
% Bottom of intergration range is eps (or realmin)
% to avoid div by 0
    sprintf('Number %g -> PiVal = %g', Number, PIval ) 
end       

function IntVal = IntFunction (w,S0,K,r,T,Model,Number,InPar)
% IntFunction is the function 'Integrated' by the PI function 
% First, set phi to name of characteristic funciton 
    if (strcmp(Model,'Gaussian'))
        phi='phiBS';
        param.sigmaBS=InPar.volBS; 
        param.mu=r-0.5*param.sigmaBS^2; %risk-neutral
        param.T=T;
    elseif (strcmp(Model,'Merton'))
        phi='phiMerton';
        param.sigmaBS=InPar.volBS; sigmaBS=param.sigmaBS;
        param.muJ=InPar.muJ;       muJ=param.muJ;
        param.sigmaJ=InPar.sigmaJ; sigmaJ=param.sigmaJ;
        param.lambda=InPar.lambda; lambda=param.lambda;
% risk-neutral drift
        param.mu=r-0.5*sigmaBS^2 ...
                    -lambda*(exp(muJ+0.5*sigmaJ^2)-1); 
        param.T=T;
    else
        disp('Improper Model')
    end
    k=log(K/S0); % dimensionless moneyness
% Similar to Black-Scholes, PI1 calculated with stock as
% numeraire and PI2 calcuted with riskless bond as numeraire
    if (Number == 1)
        IntVal=real(exp(-i*w*k).*feval(phi,w-i,param)...
                               ./(i.*w.*feval(phi,-i,param)));
    else % (Number == 2)
        IntVal=real(exp(-i*w*k).*feval(phi,w,param)...
                               ./(i*w));  
    end           
end


function c = MertonSeriesCall(K,S,T,volBS,...
                        r, d, muJ, sigmaJ, lambda);
% MertonSeriesCall based on Series Approximation where
% n is the probability of n jumps in one time period
% lambda is the intensity parameter = mean number of jumps
% in one time period 
%k=exp(muJ+0.5*sigmaJ^2)-1 = mean relative asset jump size
    N=5; % Max Number of Possible jumps in one time period
    lambdaBar=lambda*exp(muJ+0.5*sigmaJ^2);
    c=0;
% Can write following more compactly using 
% exp(muJ+0.5*sigmaJ^2)= k+1
for n=0:N
    volM = sqrt(volBS^2+n*sigmaJ/T);
    rn = r - lambda*(exp(muJ+0.5*sigmaJ^2)-1)...
        + n*(muJ+0.5*sigmaJ^2)/T;
    c = c + ((lambdaBar*T)^n / factorial(n))...
        * BlackScholesCall (K,S,T,volM,rn,0);
end
    c = c*exp(-lambdaBar*T);
end


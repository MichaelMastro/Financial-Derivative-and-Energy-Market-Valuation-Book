function LikeSum = EKFbs(parameter)
% EKF bs : Extended Kalman Filter for Black-Scholes model 
% Known Inputs: Stock Index Price, Risk Free Rate
% Observation: 1 month Call Option with (ATM) Strike = S 
% Hidden Latent Variable: volatility

% EKF appropriate for BS model as analytical equation (vega)
% is available for derivative of call price with respect 
% to volatility, i.e., linearization of the observation
% model

% Hidden volatility state follow a regressive process

global obsS obsC obsRF TimeDelta 
% Retrieve Data from Filter for Plots
global xVec hxVec MRvec volSvec

    Mu = parameter(1);	
    R = parameter(2); %Measurement Noise Covariance
    Sig = parameter(3);
    Alpha = parameter(4); 
    Beta = parameter(5); %Beta=M= State Transition Matrix
       
% Kalman filter equations are described by 
% x = volatility = C+M*S(ti-1)+B*u+w  w~N(0,Q)
% y = Call Option Price  = h(x) + v   v~N(0,R)
% h(x) = Non-Linear observation = Black-Scholes Model
% H = Gradient of BS = dc/dvol = Vega
    Strike = obsS(1);
    d1 = (log(obsS(1)/Strike)+(obsRF(1)+Mu.^2/2)*TimeDelta)...
      / (Mu*sqrt(TimeDelta));
    hxVec(1) = BlackScholesCall (Strike,obsS(1),...
                TimeDelta,sqrt(abs(Mu)),obsRF(1),0); 
    MRvec(1) = obsC(1) - hxVec(1);
    H = obsS(1)*sqrt(TimeDelta)*myNormPDF(d1); % Vega
    M = Alpha;
    C = Mu;  
    B = Beta; % Control Input 'Matrix'  
    Q = Sig*Sig*(TimeDelta);
    Pxx=Q; % Initialize State Process Covariance
    
	nsamples=length(obsC);
	LikeSum=0;%0.5*nsamples*log(2*pi);
    x=Mu;  xVec(1)=x;  varSvec(1) = 0;  MR = 0;
    
% Filter every time point
   	for i=2:nsamples 
        volS= abs((obsS(i)-obsS(i-1))...
            /(sqrt(TimeDelta)*obsS(i-1))); %Control Input = u
% Regressive State Model with volS as control input
        x = M*x+B*volS+C  ;   %x(k|k-1) State Prediction
        Pxx = M*Pxx*M'+ Q ;  %P(k|k-1) Covariance Prediction
        
        Strike=obsS(i); % at the money call option
        d1=(log(obsS(i)/Strike)+(obsRF(i)+x.^2/2)*TimeDelta)...
                / (x*sqrt(TimeDelta));
        H=obsS(1)*sqrt(TimeDelta)*myNormPDF(d1); % Vega
        hx = BlackScholesCall (Strike,obsS(i),...
            TimeDelta,x,obsRF(i),0);         %h(x)
        MR = obsC(i) - hx ;  %Measurement (Innovation)  Residual
        Pyy = H*Pxx*H' +R; %= Residual (Innovation) Covariance
        Pxy = Pxx*H';
        K = Pxy * inv(Pyy); %Kalman Gain
        K = K/4; % Damp if EKF of BS Performs Poorly when vol ~ 0
        x = x + K*(MR) ;
        % x(k|k) State Correction based on observation res.
        % Joseph Form of Covariance Correction
        Pxx = (1-K*H)*Pxx*(1-K*H)' + K*R*K' ;
        % Pxx = Pxx - K*Pyy*K';  %Pxx(k|k) Covariance Correction
        % the max. likelihood was multiplied by -1 
        % to give a minimum for fminsearch
      	LikeSum=LikeSum+0.5*log(det(Pyy))+0.5*MR'*inv(Pyy)*MR;
        % Save Information as Vector to Plot
        xVec(i)=x; hxVec(i)=hx; MRvec(i)=MR; volSvec(i) = volS;
    end
 
    % Force filter away from unreasonable values
    minx=min(xVec);
    if (Sig <=0), LikeSum=LikeSum+1e2; end
    if (R<1e-7), LikeSum=LikeSum+1e2; end
    if (minx<1e-6), LikeSum=LikeSum+1e10; end
    if (Alpha < 0), LikeSum=LikeSum+1e2; end
    if (Beta < 0), LikeSum=LikeSum+1e3; end
    if (Alpha + Beta >=1), LikeSum=LikeSum+1e2; end        
end
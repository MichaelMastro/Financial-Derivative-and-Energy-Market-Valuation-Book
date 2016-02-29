function LikeSum = UTbs(parameter)
% UTbs : Unscented Transform for Black-Scholes model 
% Known Inputs: Stock Index Price, Risk Free Rate
% Observation: 1 month Call Option with (ATM) Strike = S 
% Hidden Latent Variable: volatility
% Volatility prediction based on volatility of stock 
% (as a control input) and weighted previous volatility

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

    MRvec(1) = 0; 
    M = Alpha;
    C = Mu;  % Long Term Volatility
    B = Beta; % Control Input 'Matrix'  
    Q = Sig*Sig*(TimeDelta);
    Pxx=Q; % Initialize State Process Covariance
    
	nSteps=length(obsC(1,:));
	LikeSum=0;%0.5*nSteps*log(2*pi);
    x=Mu;  xVec(1)=x;  volSvec(1) = 0;  MR = 0;
    
    nSamp=length(x);
% large kappa will create negative sigma (volatility) points,
% which is not appropriate
    kappa=1;%3-nHid; 
    AlphaWeight = 1;
    BetaWeight= 0; % 2 = Gaussian Assumption in SUT
    lambda=AlphaWeight^2*(nSamp+kappa)-nSamp;
    Wm(1)=lambda/(nSamp+lambda); % n=1 -> W0 =2/3;
    Wc(1)=Wm(1) + (1-AlphaWeight^2 + BetaWeight);
    %Pa = blkdiag(Pxx ,Q,R);
    Wm(2:2*nSamp+1)=1/(2*(nSamp+lambda)); 
    Wc(2:2*nSamp+1)=Wm(2); % nSamp=1 -> W1 = W2 =1/6;
    
% Filter every time point %%%%%%%%%%%%%%
   	for i=2:nSteps 
% Predict
        X(1) = x';
        SigmaShift=(chol((nSamp+lambda)*Pxx))';
        X(2:nSamp+1)=X(1)+SigmaShift;  
        X(nSamp+2:2*nSamp+1)=X(1)-SigmaShift; %3:3
       
        volS= abs((obsS(i)-obsS(i-1))...
            /(sqrt(TimeDelta)*obsS(i-1))); %Control Input = u
    % Regressive State Model with volS as control input
    % Input x(k|k-1) + Sigma Points through State Prediction
        X = M*X+B*volS+C ;    
        x=sum(Wm.*X);
        %P(k|k-1) Covariance Prediction
        Pxx=(X-x)*diag(Wc)*(X-x)' + Q; % Additive Gaussian Noise Q

% Update 
        X(1) = x';
        SigmaShift=(chol((nSamp+lambda)*Pxx))';
        X(2:nSamp+1)=X(1)+SigmaShift;  
        X(nSamp+2:2*nSamp+1)=X(1)-SigmaShift; %3:3
        Strike=obsS(i);
        for j=1:(2*nSamp+1)
            Y(j) = BlackScholesCall (Strike,obsS(i),TimeDelta,...
                X(j),obsRF(i),0);  
        end
        y = sum(Wm.*Y);
   
        Pyy = (Y-y)*diag(Wc)*(Y-y)'+R; % Residual (Innovation) Cov
        Pxy = (X-x)*diag(Wc)*(Y-y)'; % Pxx*H';
        K = Pxy * inv(Pyy); %Kalman Gain
     
        MR = obsC(i) - y ;  %Measurement (Innovation)  Residual
        x = x + K*(MR) ;
        %x(k|k) State Correction based on observation res.
        %Joseph Form of Covariance Correction
        Pxx = Pxx - K*Pyy*K';  %Pxx(k|k) Covariance Correction
        %the max. likelihood was multiplied by -1 
        %to give a minimum for fminsearch
      	LikeSum=LikeSum+0.5*log(det(Pyy))+0.5*MR'*inv(Pyy)*MR;
        % Save Information as Vector to Plot
        xVec(i)=x; hxVec(i)=y; MRvec(i)=MR; volSvec(i) = volS;
    end

% Filter unstable when volatility or interest rates near zero
    minx=min(xVec); 
    if (Sig <=0), LikeSum=LikeSum+1e7; end
    if (R<1e-7), LikeSum=LikeSum+1e2; end
    if (minx<1e-6), LikeSum=LikeSum+1e12; end
    if (Alpha < 0), LikeSum=LikeSum+1e7; end
    if (Beta < 0), LikeSum=LikeSum+1e13; end
    if (Alpha + Beta >=1), LikeSum=LikeSum+1e8; end        
end
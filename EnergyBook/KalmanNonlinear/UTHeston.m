function LikeSum = UTHeston(parameter)
% UTHeston : Unscented Transform for Heston model 
% Known Inputs: Risk Free Rate
% Observation: 
        % Stock Index Price, 
        % 1 month Call Option with (ATM) Strike = S 
% Hidden Latent Variable: variance

% Call option observation model is based on 
% Fractional FFT calculation of call option via
% Heston characteristic function
% Fractional FFT is necessary to create a grid of
% call values sufficient to distinguish between sigma point
% variance values. In other words, closely spaced sigma 
% points will create call options with same value.
% This is yield unstable values for Pyy

% State Transition and Observation Model described in 
% A. Javaheri, D. Lautier, A. Galli, Filtering in Finance
    % Wilmott Magazine (2003)
% FFT of Heston for Call Option for UT described in
% J. Li, An Unscented Kalman Smoother for Volatility
% Extraction: Evidence from Stock Prices and Options,
% Computational Statistics and Data Analysis (2012)

global obsS obsC obsRF TimeDelta 
% Retrieve Data from Filter for Plots
global xVec SVec MRvecS MRvecC CVec

	PiW = parameter(1);
    Sig = parameter(2); 
    Kappa = parameter(3); 
    Theta = parameter(4);
    Rho = parameter (5);
      
% State Transition equations are described by 
% x = volatility = x + (kappa*theta - rho*sigma*(r+PiW*V)
%    -(kappa - 0.5*rho*sigma)*x)*dt + rho*sigma*ln(S(t)/S(t-1))
    % w~N(0,Q) ~ sigma*sqrt(1-rho^2)*sqrt(x*dt)
    
% y = [Index; Call Option Price]  = h(x) + v   v~N(0,R)
 % Assume Market price of Diffusive and Volatility risk are equal
    Pxx=1e-2; % Initialize State Process Covariance  
	nSteps=length(obsS(:,1));
	LikeSum=0;%0.5*nSteps*log(2*pi);
    x=xVec(1); % Kappa; 
    SVec(1)=500; CVec(1)=100; MRvecS(1:2) = 0;MRvecC(1:2) = 0;
   
    nHid=length(x);
    nSigPts=2*nHid+1;
    %Shrink Kappa Weight to avoid sigma point less than zero
    kappaWeigth=1;%3-nHid; 
    AlphaWeight = 1;
    BetaWeight= 0; % 2 = Gaussian Assumption in SUT
    lambda=AlphaWeight^2*(nHid+kappaWeigth)-nHid;
    Wm(1)=lambda/(nHid+lambda); % n=1 -> W0 =2/3;
    Wc(1)=Wm(1) + (1-AlphaWeight^2 + BetaWeight);
    %Pa = blkdiag(Pxx ,Q,R);
     Wm(2:2*nHid+1)=1/(2*(nHid+lambda)); 
     Wc(2:2*nHid+1)=Wm(2); % nHid=1 -> W1 = W2 =1/6;
     
     %Q = Sig*Sig*(TimeDelta);
     y=CVec(1)*ones(2,1);
     R=diag([1e-5 1e-5]);
     lnStm1=log(obsS(1));
     logReturn=0;
% Filter every time point %%%%%%%%%%%%%%
   	for i=2:nSteps 
% Predict      
        xPrev=x;  
        X(1) = x';
        SigmaShift=(chol((nHid+lambda)*real(Pxx)))';
        X(2:nHid+1)=X(1)+SigmaShift; %2:2
        X(nHid+2:nSigPts)=X(1)-SigmaShift; %3:3
        
    % Input Sigma Points into State Prediction   
        X = X+(Kappa*Theta-Rho*Sig*(obsRF(i)+PiW*X)...
            -(Kappa-0.5*Rho*Sig)*X)*TimeDelta+ Rho*Sig*logReturn;
        logReturn=log(obsS(i)/obsS(i-1));
        x=sum(Wm.*X); 
        %P(k|k-1) Covariance Prediction
        Q = (Sig*sqrt(1-Rho^2)*sqrt(x*TimeDelta))^2 + 1e-9;     
        Pxx=((X-x)*diag(Wc)*(X-x)') + Q; % Additive Gaussian Noise Q
        
 % Update 
        X(1) = x';
        SigmaShift=(chol((nHid+lambda)*real(Pxx)))';
        X(2:nHid+1)=X(1)+SigmaShift; %2:2
        X(nHid+2:nSigPts)=X(1)-SigmaShift; %3:3
                
        Y(1,:) = lnStm1 + (obsRF(i)+PiW.*X(:) - 0.5*X(:))*TimeDelta;
    % Call option via Fractional FFT of Heston Characteristic func.
   
        for j=1:nSigPts
           Y(2,j) = CalcHestonCall (parameter, obsS(i-1),...
                    obsS(i-1), obsRF(i), X(j), TimeDelta);     
        end
        lnStm1=log(obsS(i));
 
        y = sum((repmat(Wm,2,1).*Y),2);
        deltaY=Y-repmat(y,1,nSigPts) ;  


        R =  (sqrt(x*TimeDelta))^2+ 1e-9;
        Pyy = (deltaY*diag(Wc)*deltaY')+R; % Residual (Innovation) Cov
        Pxy = ((X-x)*diag(Wc)*deltaY'); % Pxx*H';
        K = Pxy * inv(Pyy); %Kalman Gain
     
        MR = [log(obsS(i)); obsC(i)] -y;  %Measurement (Innovation)  Residual
        x = x + K*(MR);% log(obsS(i));
        
        %x(k|k) State Correction based on observation res.
        %Joseph Form of Covariance Correction
        Pxx = Pxx - K*Pyy*K';  %Pxx(k|k) Covariance Correction
  
        %the max. likelihood was multiplied by -1 
        %to give a minimum for fminsearch
      	LikeSum=LikeSum+0.5*log(det(Pyy))+0.5*MR'*inv(Pyy)*MR;
        % Save Information as Vector to Plot
        xVec(i)=x;  
        MRvecS(i) = MR(1); MRvecC(i) = MR(2); 
        SVec(i) = y(1); CVec(i) = y(2);
    end 
end







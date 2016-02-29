% OUloglikelihoodfn returns negative of log Maximum Likelihood 
% Estimationof process parameters (mu, lambda, sig, R) fit 
% to observation data with the Kalman filter
function LikeSum = OUloglikelihoodfn(parameter)
    global Obs; %use of global avoids passing large arrays
    global TimeDelta; 
    global xVec; % Use to Save Final Predicted Price Series
% Only mu, lambda, sig, and R variables are optimized by 
% fminsearch in calling function based on likelihood metric
% that is output of OUloglikelihoodfn  
	mu=parameter(1);
	lambda=parameter(2);
	sig=parameter(3);
    R=parameter(4); %Measurement Noise Covariance
    
%Kalman filter equations are described by 
%x=S(ti)=C+M*S(ti-1)+w  w~N(0,Q)
%y=Obs  =Hx+v           v~N(0,R)
    H=1; 
  %Directly Observing Price -> Observation Matrix = 1 (or Identity)
    M=exp(-lambda*TimeDelta);
    C=mu.*(1-M); %=mu.*(1-exp(-lambda*TimeDelta));
    %Q=State Space Process Variance
%if (lambda==0) %include if-else if searching for process w/ lambda~0
        Q= sig*sig*( (1-exp(-2.*lambda.*TimeDelta)) ./ (2.*lambda));
    %else
    %   Q= sig*sig*(TimeDelta);
%end
    
	nsamples=length(Obs);%size(Obs,2);
	LikeSum=0.5*nsamples*log(2*pi);
    x=inv(H)*Obs(1);
    xVec(1)=x;
    % Initialize Covariance to Unconditional Variance of State Process
    P=inv(H)*Q*inv(H'); 
% approach avoids saving every a-priori and a-posteriori variable
   	for i=2:nsamples  
        x = M*x+C;       %x(k|k-1) State Prediction
        P = M*P*M'+ Q;   %P(k|k-1) Covariance Prediction
   
        MR = Obs(i) - H*x ; %Measurement (Innovation)  Residual
		V = H*P*H' + R;    %Residual (Innovation) Covariance
        K = P*H' * inv(V); %Kalman Gain
        
        x = x + K*(MR); 
        %x(k|k) State Correction bassed on observation res.
        P = P - K*H*P;  %P(k|k) Covariance Correction
        %the max. likelihood was multiplied by -1 
        %to give a minimum for fminsearch
      	LikeSum=LikeSum+0.5*log(det(V))+0.5*MR'*inv(V)*MR;%+1e-10;
        xVec(i)=x;
    end
end
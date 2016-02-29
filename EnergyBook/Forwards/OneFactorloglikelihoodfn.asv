% OneFactorloglikelihoodfn returns negative of 
% log Maximum Likelihood Estimation of process parameters 
% (alpha-real-world, alpha-risk free lambda, sig) fit to 
% observation data with the Kalman filter
function LikeSum = OneFactorloglikelihoodfn(parameter)
    global Obs; %use of global avoids passing large arrays
    global TimeDelta; 
    global TimeMat;
    global LnS;
    global MeasErr;
    
	alpha=parameter(1); kappa=parameter(2); 
    sig=parameter(3); alphaStar=parameter(4); 
    Rvar=parameter(5:8); R=diag(Rvar);%Measurement Noise Covariance

    Ncontracts=length(TimeMat);
%Kalman filter equations are described by 
%x=ln(S(ti))=C+M*(ln(S(ti-1))+w  w~N(0,Q)
%y=Obs  =D+Hx+v           v~N(0,R)
    kappa=max(0.01,kappa); %Avoid Neg Reversion
    %relate spot price to forwards at different maturities  
    D=(alphaStar.*(1-exp(-kappa.*TimeMat))...
        +(sig*sig./(4*kappa))*(1-exp(-2*kappa*TimeMat)))';
    H=exp(-kappa*TimeMat)'; %4x1 Matrix
    
    M=exp(-kappa*TimeDelta); %M=1x1 Matrix in 1 factor Model
    C=alpha*(1-M); %=alpha.*(1-exp(-kappa*TimeDelta));
    Q= sig*sig*( (1-exp(-2.*kappa.*TimeDelta)) ./ (2.*kappa));
      
	nsamples=length(Obs(:,1));
    nseries=length(Obs(1,:));
	LikeSum=0.5*nsamples*log(2*pi);
    
    x=Obs(1,1)/H(1);%inv(H)*Obs(1,:);
 %Initialize Covariance to Unconditional Variance of State Process
    P=Q;
   	for i=2:nsamples  
        x = M*x+C;       %x(k|k-1) State Prediction x=4x1 Matrix
        P = M*P*M'+ Q;   %P(k|k-1) Covariance Prediction

        MR = Obs(i,:)' - H*x-D;  %Measurement (Innovation) Residual
		V = H'*P*H + R;    %Residual (Innovation) Covariance
        K = P*H' * inv(V); %Kalman Gain 
        
        x = x + K*(MR) ;
        %x(k|k) State Correction bassed on observation res.
        P = P - K*H*P ; %P(k|k) Covariance Correction
        %the max. likelihood was multiplied by -1 
        %to give a minimum for fminsearch
      	LikeSum=LikeSum+0.5*log(det(V))+0.5*MR'*inv(V)*MR;%+1e-10; 
        LnS(i)=x;
        MeasErr(i,:)=MR;
    end
    for i = 1:8 % brute force technique to avoid negative terms
        if (parameter(i)<=0)
            LikeSum=LikeSum+1e5;
        end
    end
end



function LikeSum = A3Mloglikelihoodfn(param)%,knownParam) 
% A3Mloglikelihoodfn returns negative of 
% log Maximum Likelihood Estimation of 3 Component Affine Model
% via Kalman Filter. Jump were pre-filtered in calling program

global tau
global Obs
global knownParam
global latent
global MRvec

kappaY = knownParam(1);
lambda = knownParam(2);
muJ = knownParam(3);
sigJ = knownParam(4);
dt = knownParam(5);
nObs=knownParam(6);
nsamples=knownParam(7);

muChi = param(1);
kappaChi = param(2);
sigChi = param(3);
muE = param(4);
kappaE = param(5);
sigE = param(6);
rho  = param(7); 
phiChi =param(8);
Rvar=param(9:13);   % Variance terms on diagonal of
R=diag(Rvar);           % Measurement Noise Covariance

% Kalman filter equations are described by 
% x(t) = M*x(t-1) + C + w   w~N(0,Q)
% z = Obs = Hx + d + v    v~N(0,R)

    H = ones(nObs,2); % H + d = Observation Model
 %  H(:,1) = exp(-kappaChi*tau(1,:)); %
% Measurement Noise Covariance is small but non-zero
    %R = diag(0.1*ones(nObs,1)); 
% Nomikos suggests Zero Measurement Noise for Spot
% but here we will simply make smaller % R(1,1)=0.001;
  
% State Transition Matrix M
% Chi (short term variation) exponential decay (exp(-kappaChi*dt))
% Epsilon (long term rate)  drifts at rate (muE-0.5*sigE^2*dt)
% Linear M + Constant C = State Transition Model
    M = [exp(-kappaChi*dt), 0; 0, 1];
    C = [0; muE-0.5*sigE^2*dt];
    % Q = State Space Process Variance
    % include if-else if searching for process w/ kappaChi==0
    if (kappaChi==0) 
        Q=sig*sig*(dt);
        Q(1,1)=sigChi*sigChi*dt;  
        Q(1,2)=rho*sigE*sigChi*dt;
        Q(2,1)=Q(1,2);
        Q(2,2)=sigE*sigE*dt;
    else
        Q(1,1)=sigChi*sigChi*( (1-exp(-2.*kappaChi.*dt))...
            ./ (2.*kappaChi));  
        Q(1,2)=rho*sigE*sigChi*( (1-exp(-kappaChi.*dt))...
            ./ kappaChi);
        Q(2,1)=Q(1,2);
        Q(2,2)=sigE*sigE*dt;
    end
    
	LikeSum=0.5*nsamples*log(2*pi);
    % Create Crude Initial Estimation of Latent Variable
    x=[Obs(1,1)-Obs(end,1) ;Obs(end,1)]; 
%Initialize Covariance to Unconditional Variance of State Process
    P=Q;
    
%%%% Step Forward in Time %%%%%%%%%%%%%%%%%%%%%
%approach avoids saving every a-priori and a-posteriori variable
for i=2:nsamples  
    H(:,1) = exp(-kappaChi*tau(i,:));
    % Second column of Observation Model (that multiples Epsilon)
    % is always equal to unity % H(:,2) = ones(nObs,1); 

% Integrate Jump Component by summing Rectangles
    djump=zeros(1,nObs);
    maxj=5; jcount=linspace(0,maxj,maxj);
    deltaTau=tau(i,:)/maxj;
    for j=1:maxj
        % tau mid-point of jth rectangle
        TauChop=(j-0.5)*deltaTau; 
        % jth Heigth at mid-point
        djump=djump+(exp(muJ*exp(-kappaY*TauChop)...
            +0.5*sigJ^2*exp(-2*kappaY*TauChop))-1);
    end
    djump=lambda*djump.*deltaTau; % lambda = jump frequency
 % (Linear) H + (Constant) d = Observation Model maps hidden
 % (latent) variables X, eps into observation space
    d=(sigChi^2/(4*kappaChi))*(1-exp(-2.*kappaChi.*tau(i,:)))...
        -((phiChi-rho*sigE*sigChi)/kappaChi)...
           *(1-exp(-kappaChi.*tau(i,:))) + muE*tau(i,:) + djump;
% Predict
        x = M*x+C;      %x(k|k-1) State Prediction
        P = M*P*M'+ Q;  %P(k|k-1) Covariance Prediction
%Update
        MR = Obs(:,i) - H*x-d';%Measurement (Innovation) Residual
		V = H*P*H' + R;    %Residual (Innovation) Covariance
        K = P*H' * inv(V); %Kalman Gain
        %x(k|k) State Correction based on observation Residual    
        x = x + K*(MR);       
        P = P - K*H*P;  %P(k|k) Covariance Correction
  
% negative of maximum likelihood gives minimum for fminsearch  
        LikeSum=LikeSum+0.5*log(det(V))+0.5*MR'*(V\MR); 
        MRvec(i,:) = MR;  % store globally to plot later    
        latent(i,:)=x; % store globally to plot later
end
% Crude technique to prevent unreasonable parameters
if ((sigE<=0)||(sigChi<=0)||(kappaE<0)||(kappaChi<0))
       LikeSum=LikeSum+1e7;
end 
   for i = 9:13 % brute force technique to avoid negative variance
        if (param(i)<=0)
            LikeSum=LikeSum+1e7;
        end
    end
end
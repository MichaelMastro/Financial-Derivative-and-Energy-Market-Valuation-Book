% CYloglikelihoodfn 2-factor model returns negative of 
% log Maximum Likelihood Estimation of process parameters 
% fit to (futures contract) data with the Kalman filter
% dS=(mu-CY)Sdt+sig1*S*dz1
% dCY=kappa(alpha-CY)dt+sig2*dz2
function LikeSum = CYloglikelihoodfn(parameter)
% Four Filtering Techniques are available: flag =
% 0='Normal-Joseph';%1='scalar';%2='Schmidt-SR';%3='Carraro-S/SR';

% 0='Normal-Joseph': uses More Symmetric Joseph Covariance Correction

% 1='scalar': when measurents cov terms = 0 -> Sequential scaler 
% Filtering allows division instead of inversion 

% 2='Schmidt-SR' only used if zero correlation between measurements.
% Improves inversion by Doubling the Precision 
% of ill-conditioned Innovation Covariance Matrix
% Highly correlated measurements precisely measured, i.e.,
% with very small variance terms on the diagonal of the Measurement 
% Noise Covariance matrix, create a situation where all terms 
% in the Innovation Covariance Matrix are similar. This leads
% to ill-conditioned matrices during the inversion step.

% 3='Carraro-S/SR'
% Propagates the squared covariance matrix itself through the filter, 
% thus losing the increase in numerical precision
% Forces the calculation of a positive semi-definite covariance matrix

    global Obs ; % use of global avoids passing large arrays
    global TimeDelta; % generally global variables are bad programming
    global TimeMat;
    global LnS CY MeasErr flag;

    r=0.06; % Constant interest Rates -> Futures=Forward Contract
	alpha=parameter(1); kappa=parameter(2); 
    sig1=parameter(3); lambda=parameter(4);
    sig2=parameter(5); rho=parameter(6);
    mu=parameter(7);
    Rvar=parameter(8:11);   % Variance terms on diagonal of
    R=diag(Rvar);           % Measurement Noise Covariance
   
    Ncontracts=length(TimeMat);
    % lnF=lnS + CY*(1-exp(-kappa*dt))/kappa + D
    % Measurement Equation
    % y[4x1] = H[4x4] * x[2x1] + D[4x1]
    % relate spot price to forwards at different maturities 
    H=[ones(Ncontracts,1),(-(1-exp(-kappa*TimeMat))/kappa)];
    D=((r+0.5*sig2^2/kappa^2-alpha+sig2*lambda/kappa...
                -sig1*sig2*rho/kappa)*TimeMat)...
        +(0.25*sig2^2*(1-exp(-2*kappa*TimeMat))/kappa^3)...
    +((alpha-sig2*lambda/kappa+sig1*sig2*rho/kappa-sig2^2/kappa^2)...
            *(1-exp(-kappa*TimeMat)))/kappa;

    % lnS(t)=lnS(t-1)-CYdt+(mu-0.5*sig1*sig1)dt
    % CY(t)=CY(t-1)*exp(-kappa*dt)+alpha(1-exp(-kappa*dt)
    % Transition Equation    
    % x[2x1] = M[2x2] * x[2X1] + C[2x1] + w  w~N(0,Q)
    M=[1, -TimeDelta ; 0, exp(-kappa*TimeDelta) ]; %M=2x2 Matrix 
    C=[TimeDelta*(mu-0.5*sig1^2); alpha*(1-exp(-kappa*TimeDelta))];
    
    Qv1= sig1*sig1*TimeDelta;
    Qv2= sig2*sig2*( (1-exp(-2.*kappa.*TimeDelta)) ./ (2.*kappa));
    Qc= rho*sig1*sig2*( (1-exp(-kappa.*TimeDelta)) ./ kappa);
    Q=[Qv1, Qc; Qc, Qv2];
    
	nsamples=length(Obs(:,1));
    nseries=length(Obs(1,:));
	LikeSum=0.5*nsamples*log(2*pi);
    
    x=[Obs(1,1)/H(1,1); 0];%inv(H)*Obs(1,:);
    %Initialize Covariance to Unconditional Variance of State Process
    P=[Qv1 0; 0 Qv2]; 
    I2=eye(2); I4=eye(4);
    
switch flag
case {'1','scalar'}
   % flag='scalar';
nmeasurements=length(Obs(1,:));
for i=2:nsamples 
    % update state and covariance predictions as a vector
        x = M*x+C;       %x(k|k-1) State Prediction x=4x1 Matrix
        P = M*P*M'+ Q;   %P(k|k-1) Covariance Prediction
    % Individually filter each futures contract    
    for j=1:nmeasurements %step through 4 measurements / time step 
        zpred=H(j,:)*x-D(j); %scalar
        MR(j)=Obs(i,j)-zpred; %scalar
        HP=H(j,:)*P; %HP[1x2]=H(j)[1x2] P[2x2]
        V(j)=[[HP]*H(j,:)'+R(j,j)]; %[1x4]
        invV(j)=1/V(j); %scalar-> do not need invert matrix
        K=(HP)'*invV(j);
        x=x+K*MR(j);
        P=P-K*(HP);
    end
    % put into form for maximum likelihood calculation
    dInvV=diag(invV);
    % DiagonalV=diag(V); DetDiagonalV=det(DiagonalV)
    ProdDiag=prod(V);
    % determinant of diagonal matrix is product elements
    LikeSum=LikeSum+0.5*log(ProdDiag)+0.5*MR*dInvV*MR';
    LnS(i)=x(1);
    CY(i)=x(2);
    MeasErr(i,:)=MR;
end
  
case {'2','S-SR','Schmidt-SR'}
    % flag='Schmidt-SR';
    S=(chol(P))'; %S(0|0) 
    Q1=(chol(Q))';
    R1=(chol(R))';
    for i=2:nsamples 
        x = M*x+C ; %x(k|k-1)=M*x(k-1|k-1)
        A1=[S'*M';Q1']'; %[2x4] matrix:
        % A1 is rectangular and not Triangular
        % convert A1 into upper triangle S such that 
        % A1*A1'=S*T*T'*S'= S*S'= M*P*M'+ Q 
        [T S]=qr(A1',0); %T*S=A1' and S=T'*A1'
        %Matlab function qr factors A1 into upper triangular matrix S
% and orthogonal matrix T that is discarded
% qr could be coded with a givens or householder function
        S=S'; %lower triangle
        %A1*A1'=S*S'= M*P*M'+ Q 
        %triangular S can be propagated through filter
        A2=[S'*H';R1']'; %(k|k-1)
        V=A2*A2';
        B=chol(V);%upper triangle% V=B'*B
        invB=inv(B);
        invV=invB*invB';%=inv(B)*inv(B')%=inv(V)
        detV=det(B)^2;%det(V)=det(B)*det(B')=det(B)*det(B)     
        
        K=S*S'*H'*invB*invB'; %Kalman Gain 
        MR = Obs(i,:)' - H*x-D;  %Innovation Residual
        %x(k|k) State Correction based on observation res.
        x = x + K*(MR) ;
        gamma=inv(I4+invV*R); %  %scalar 
        S=(I2-K*gamma*H)*S;%S(k|k) %scalar 
        %Never necessary to calculate Covariance P =S*S' 
        LikeSum=LikeSum+0.5*log(detV)+0.5*MR'*invV*MR;
        LnS(i)=x(1);
        CY(i)=x(2);
        MeasErr(i,:)=MR;
    end


case {'3','Carraro-S/SR','S/SR'}
    % flag='Carraro-S/SR';
    % Carraro Squaring/Square root Kalman filter
% The product of A(j)*A(j)'=P j=1,2,3 should always be pos semi-def 
    Q1=(chol(Q))';
    R1=(chol(R))';
    %size(R1)
    for i=2:nsamples 
        S=(chol(P))'; %S(0|0)  or %S(k|k) 
        x = M*x+C ; %x(k|k-1)=M*x(k-1|k-1)+C
        
        A1=[S*M';Q1']';
        P=A1*A1' ;%P(k|k-1)
        S=(chol(P))'; %S(k|k-1)
        A2=[S'*H';R1']';
        V=A2*A2'; %%Residual (Innovation) Covariance
        invV=inv(V);
        K = P*H' * invV ;%P*H' * inv(V) %Kalman Gain 
        MR = Obs(i,:)' - H*x-D;  %Measurement (Innovation)  Residual
        % x(k|k) State Correction based on observation res.
        x = x + K*(MR) ;
        A3=[S'*(I2-K*H); R1'*K']';
        P=A3*A3'; %P(k|k)
        LikeSum=LikeSum+0.5*log(det(V))+0.5*MR'*invV*MR;%+1e-10;
        LnS(i)=x(1);
        CY(i)=x(2);
        MeasErr(i,:)=MR;
    end
    
otherwise % default to {'0','Normal-Joseph'}
   % flag='Normal-Joseph';
    for i=2:nsamples 
        x = M*x+C  ;     % x(k|k-1) State Prediction x=4x1 Matrix
        P = M*P*M'+ Q ;  % P(k|k-1) Covariance Prediction
        MR = Obs(i,:)' - H*x-D ; %Measurement (Innovation)  Residual
        V = H*P*H' + R;    % Residual (Innovation) Covariance
        invV=inv(V);
        K = P*H' * invV; % P*H' * inv(V) %Kalman Gain 
        % x(k|k) State Correction based on observation res.
        x = x + K*(MR) ;
        % P(k|k)  Covariance Correction
        % -> use More Symmetric Joseph Covariance Correction
        P=(I2-K*H)*P*(I2-K*H)'+K*R*K';
        LikeSum=LikeSum+0.5*log(det(V))+0.5*MR'*invV*MR;%+1e-10;
        LnS(i)=x(1);
        CY(i)=x(2);
        MeasErr(i,:)=MR;
    end
end
    for i = 1:11 % brute force technique to avoid negative terms
        if (parameter(i)<=0)
            LikeSum=LikeSum+1e5;
        end
    end
end
%-------------------------------------------------------------



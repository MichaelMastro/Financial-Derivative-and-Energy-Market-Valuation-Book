function LikeSum = PFHeston(parameter)
% PFHeston : Particle Filter for Heston model 
% FFT calculation slows down entire calculation
% so limit number of particles to keep calculation time
% reasonable
% Known Inputs: Risk Free Rate
% Observation: 
        % Stock Index Price, 
        % 1 month Call Option with (ATM) Strike = S 
% Hidden Latent Variable: variance

% Call option observation model is based on 
% Fractional FFT calculation of call option via
% Heston characteristic function
% Fractional FFT is necessary to create a grid of
% call values sufficient to distinguish between particles
% variance values. In other words, closely spaced  
% particles will create call options with same value.

% More detailed and robust particle filters are
% available from Nando de Freitas (Berkeley) and 
% Rudolph van der Merwe (OGI) 

global obsS obsC obsRF TimeDelta 
% Retrieve Data from Filter for Plots
global xVec SVec MRvecS MRvecC CVec WeightVec PartVec

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
	LikeSum=0; %0.5*nSteps*log(2*pi);

    x=xVec(1); 
    SVec(1)=500; CVec(1)=100; MRvecS(1:2) = 0;MRvecC(1:2) = 0;
   
    nHid=length(x);
     y=CVec(1)*ones(2,1);
     R=diag([1e-2 1e-2]);
     lnStm1=log(obsS(1));
     logReturn=0;
     
     Npart=11; %%% Number of Particles
     PartX = x *ones(1,Npart);
      
% Filter every time point %%%%%%%%%%%%%%
   	for i=2:nSteps 
        
        Q = (Sig*sqrt(1-Rho^2)*sqrt(x*TimeDelta))^2 + 1e-9; 
        PartX=PartX+sqrtm(Q)*randn(1,Npart);
% Predict      
        xPrev=x;        
    % Input Particles into State Prediction   
        PartX = PartX+(Kappa*Theta-Rho*Sig*(obsRF(i)+PiW*PartX)...
            -(Kappa-0.5*Rho*Sig)*PartX)*TimeDelta+ Rho*Sig*logReturn;
        PartX=abs(PartX); % Prevent Negative Variance
        logReturn=log(obsS(i)/obsS(i-1));     
% Update    
        PartY(1,:) = lnStm1 +...
            (obsRF(i)+PiW.*PartX(:) - 0.5*PartX(:))*TimeDelta;
    % Call option via Fractional FFT of Heston Characteristic func.   
        for j=1:Npart
           PartY(2,j) = CalcHestonCall (parameter, obsS(i-1),...
                    obsS(i-1), obsRF(i), PartX(j), TimeDelta);     
        end
        lnStm1=log(obsS(i));
            
        yObs = [log(obsS(i)); obsC(i)];
        deltaY=PartY-repmat(yObs,1,Npart) ; 
        R(1,1) =  (x*TimeDelta)+ 1e-9;
   
        for j=1:Npart
            w(j)=exp(-0.5*deltaY(:,j)'*inv(R)*deltaY(:,j))+1e-15;
        end
        w=w/sum(w); % Normalize Weights 
        WeightVec(i,:)=w; PartVec(i,:)=PartX; % Save for Plotting
        yMean=sum(repmat(w,2,1).*PartY,2);
        
        x=sum(w.*PartX); 
        deltaYY=PartY-repmat(yMean,1,Npart) ; 
        Pyy = (deltaYY*diag(w)*deltaYY')+R; % Residual (Innovation) Cov
        
        %Measurement (Innovation)  Residual
        MR = [log(obsS(i)); obsC(i)] -yMean;  
        %the max. likelihood was multiplied by -1 
        %to give a minimum for fminsearch
      	LikeSum=LikeSum+0.5*log(det(Pyy))+0.5*MR'*inv(Pyy)*MR;
        
        NewIndex = Resample(w);
        % [OldIndex=1:NumParticles; OldWeights=w] 
        % to [NewIndex ; NewWeights=1/NumPart.,...,1/NumPart.]
        PartX(:)=PartX(NewIndex);
        
        % Save Information as Vector to Plot
        xVec(i)=x;  
        MRvecS(i) = MR(1); MRvecC(i) = MR(2); 
        SVec(i) = yMean(1); CVec(i) = yMean(2);
    end 
end

function NewIndex = Resample(w)
% Resample resamples [OldIndex=1:NumParticles; OldWeights=w] 
% to [NewIndex ; NewWeights=1/NumPart.,...,1/NumPart.]
% e.g., NewIndex = [1 1 4 5 5 6 6 6 6]
% Excellent comparison of resampling schemes available from
% Nando de Freitas (Berkeley) and Rudolph van der Merwe (OGI) 
    L=max(size(w));
    OldIndex=1:L;
    NewPart=zeros(1,L);
    CDFunc=cumsum(w);   

    u = sort(rand(1,L)); 
    % Nando de Freitas and Rudolph van der Merwe suggest 
    % u=fliplr(cumprod(rand(1,L).^(1./(L:-1:1))));
    j=1;
    for counter=1:L
        while(u(counter) > CDFunc(j))
            j=j+1;
        end
        NewPart(j)=NewPart(j)+1;
    end
    
    ind=1;
    for counter=1:L
        if (NewPart(1,counter)>0)
            for j=ind:ind+NewPart(counter)-1
                NewIndex(j) = OldIndex(counter);
            end
        end
        ind=ind+NewPart(counter);
    end                
end





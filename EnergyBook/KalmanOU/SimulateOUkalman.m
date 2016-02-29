%SimulateOUkalman simulates a discrete mean reversion process (S)
%Noise in corresponding observation (Obs) 
%arises from errors in data feeds
function SimulateOUkalman (mu,lambda,  Dsig, Szero, steps)
    global S; global Obs; global TimeDelta;
    %Expected = Szero.*exp(-lambda.*time)+mu.*(1-exp(-lambda*time)); 	    
    M=exp(-lambda*TimeDelta); 
    C=mu.*(1-M); %=mu.*(1-exp(-lambda*TimeDelta));
	
SD= Dsig*sqrt( (1-exp(-2.*lambda.*TimeDelta)) ./ (2.*lambda));
Q=SD.^2; %covariance of unobservable process

	R=0.01; 
%assume (co)variance of observable process is small but non-zero
	H=1; %We are directly observing spot price -> H(matrix)=1
    
	noise_w=SD.*randn(1,steps);
    noise_v=sqrt(R).*randn(1,steps);
       
    S(1) = Szero;
    Obs(1) = H*S(1)+noise_w(1);
    
	for i=2:steps
        %S= Expected + SD.*randn(1,steps);
        S(i)=C+M*S(i-1)+noise_w(i);%x=S(ti)=C+M*S(ti-1)+w    w~N(0,Q)
        Obs(i)=H*S(i)+noise_v(i); %Z=Obs=Hx+v                v~N(0,R)
	end
end
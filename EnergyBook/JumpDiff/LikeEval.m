function LH = LikeEval(parameter)
% LikeEval returns negative of CHI SQUARE (lambda, a,b) fit to data
% or negative of Multinomial Maximum Likelihood

    global M1 M2 % mean and variance
    global fData Centers lengthLD dt 
    global LSflag %1 for LS; 2 for MME
% lambda, a,b are independent variable
    lambda=parameter(1); 
    a=parameter(2); 
    b=parameter(3);
% mu and sigma of dependent on (lambda, a,b) 
% via M1 (mean) and M2 (variance)
% see F. Hanson, J.J.Westman, 
% Jump-Diffusion Stock Return Models in Finance:
% Stochastic Process Density with Uniform Jump Amplitude
    muJump=(a+b)/2;
    sigJ2=(b-a)^2/12;
    muSig2=(M1-muJump*lambda*dt)/dt;
    sig=sqrt((M2-(sigJ2+muJump^2)*lambda*dt)/dt);
    mu=muSig2+0.5*sig^2;
% uniform jump-diffusion log return probability density 
% Equation derived by D. Synowiec, 
% Comp. Math w/ App. 56, 2120 (2008), %uses 
% a normalized PDF and CDF which is calculated 
% with myNormPDF and myNormCDF
% PDF and CDF functions are not in standard Matlab package

    fjd= ((1-lambda*dt)/(sig*sqrt(dt)))...
    *myNormPDF((Centers-(mu-0.5*sig*sig)*dt) / (sig*sqrt(dt)))...
    +(lambda*dt/(b-a))*...
        ( myNormCDF((Centers-a-(mu-0.5*sig*sig)*dt)/(sig*sqrt(dt)))...
            - myNormCDF((Centers-b-(mu-0.5*sig*sig)*dt)/(sig*sqrt(dt))) );
% Normalize so total frequency(fjd)=length(Energy or Stock Prices) 
if (LSflag==1)
    fjd=lengthLD*fjd/sum(fjd);
    % Unweighted Chi Square approach: Fairly Stable 
    % and give Reasonable Results
    LH=sum((fjd-fData).^2); 
    % Hanson and Westman suggest a weighted Chi Square approach
    % fvar=abs(lengthLD*fjd.*(1-fjd/lengthLD).^2);
    % w=(1./fvar)./sum(1./fvar);   
    % weighted Can be unstable: Try alternate Weighting scheme %w=(fjd); 
    % LH=sum(w.*(fjd-fData).^2);
else
% Multinomial Maximum Likelihood
% Derived in Floyd B. Hanson, John J. Westman and Zongwu Zhu, 
% "Maximum Multinomial Likelihood Estimation of Market Parameters
% for Stock Jump-Diffusion Models, in Mathematics of Finance"
    LH=-sum((fData.*log(fjd)));
end
end
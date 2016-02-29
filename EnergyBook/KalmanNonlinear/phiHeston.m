function [ phiVal ] = phiHeston (w,Param)
% phiMerton is characteristic function for Merton's 
% jump-diffusion model
    PiW = Param.PiW;
    Sig = Param.Sig;
    Kappa = Param.Kappa; 
    Theta = Param.Theta;
    Rho = Param.Rho; 
    r = Param.r;
    V = Param.V;
    T = Param.T;
    
    KappaStar=(Kappa+PiW)-i*w*Rho*Sig;
    psiW=0.5*(i*w+w.^2);
    Gamma=sqrt(KappaStar.^2+2*Sig^2*psiW);
    A=(Kappa.*Theta./Sig.^2)...
        *(2.*log(1- ((Gamma-KappaStar).*(1-exp(-Gamma.*T)))./(2.*Gamma))...
        + (Gamma-KappaStar).*T);
    B=2*psiW.*(1-exp(-Gamma.*T))...
        ./ (2.*Gamma-(Gamma-KappaStar).*(1-exp(-Gamma.*T)));
    phiVal=exp(i*w*r*T - A - B*V);
end
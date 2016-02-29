function C = CalcHestonCall (parameter, S0,K, r, V, T)
% CalcHestonCall calls FFToption function for fractional
% FFT calculation of Heston call price
% Setting the parameter N (to any number) tells FFToption
% to use the fractional FFT
% The unscented transform calculates option values at 
% closely spaced sigma points. The Fractional FFT option algorithm 
% provides a fine grid at a reasonable cost

    CharFunc='phiHeston';
    Param.PiW = parameter(1);
    Param.Sig = parameter(2); 
    Param.Kappa = parameter(3); 
    Param.Theta = parameter(4); %Beta=M= State Transition Matrix
    Param.Rho = parameter (5);
    Param.r= r;
    Param.V =V;
    Param.T=T;

     wEnd=500;
    alpha=1.25;
    N=2^8;  
          
    [callFFT, k ] = FFToption(CharFunc, Param, alpha, wEnd , N);   

    Kvec=S0*exp(k); %K is normalized to S0 in FFT calc

    C = interp1(Kvec,callFFT,K);
    C=C*S0; 
end
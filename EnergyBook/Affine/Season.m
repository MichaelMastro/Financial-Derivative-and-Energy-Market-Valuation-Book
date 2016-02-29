function valSeason = Season (par,t)
% Two Sinusoidal function called by SineFit and Affine3Model
    amp1=par(1);
    phase1=par(2);
    amp2=par(3);
    phase2=par(4);
    valSeason=amp1*sin(2*pi*(t+phase1))...
        +amp2*sin(4*pi*(t+phase2));   
end
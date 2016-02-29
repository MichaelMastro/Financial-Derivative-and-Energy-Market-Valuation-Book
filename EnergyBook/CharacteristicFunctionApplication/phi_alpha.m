function phi= phi_alpha(omega,p); 
% phi_alpha returns characteristic function 
% of alpha stable  distribution
    t=p.t;
    mu=p.mu;
    alpha=p.alpha;
    beta=p.beta;
    sigma=p.sigma;
    if (alpha==1)
        phi=i*mu*omega...
            -sigma.*(abs(omega))...
            .*(1+i*beta*sign(omega)*(2/pi) ...
            .*(log(abs(eps+omega))));
        %add eps to avoid log(0) 
    else %alpha not 1
        phi=i*mu*omega...
            -sigma.^alpha.*(abs(omega)).^alpha...
            .*(1-i*beta*sign(omega)*tan(pi*alpha/2));
    end

    phi=exp(t*phi);
end
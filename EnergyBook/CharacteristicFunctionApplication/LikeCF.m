function LH = LikeCF(pnew)
% LikeCF calculates negative likelihood that parameter
% set pnew creates an alpha-stable distribution that 
% matches the experimental data PDFexp
% LikeCF is called with fminsearch function within
% FitAlphaDist

    global PDFexp
    global xout
    global TimeDelta
    global PDFtheo

    p.alpha=pnew(1);
    p.beta=pnew(2);
    p.mu=pnew(3);
    p.sigma=pnew(4);
    p.t=TimeDelta;
    LH=0;
    
    %Quick Crude approach to bound parameters
    if (p.alpha>2), p.alpha=2; LH=1e9+1e9*i; end
    if (p.alpha<1), p.alpha=1; LH=1e9+1e9*i;end
    if (p.beta>1), p.beta=1;LH=1e9+1e9*i; end
    if (p.beta<-1), p.beta=-1;LH=1e9+1e9*i; end   
    if (p.sigma<0.01), p.sigma=0.01;LH=1e9+1e9*i; end    
    
    OmegaEnd=round(2/(p.sigma*p.t*(p.alpha^2)));
    omega=linspace(0,OmegaEnd,1000)';
    DeltaOmega=omega(2)-omega(1);
    phin=phi_alpha(omega,p);
    
    Nx=length (xout);
    PDFtheo=zeros(Nx,1);
    
    for k=1:Nx
        y=exp(-i*xout(k)*omega).*phin;
        y=real(y);
        I=sum(y)-0.5*(y(1)+y(end)); % trapezoidal rule
        PDFtheo(k)=I*DeltaOmega/pi; %output
    end
    PDFtheo=PDFtheo/sum(PDFtheo); 
    LH=LH+sum((PDFexp(:)-PDFtheo(:)).^2);
end
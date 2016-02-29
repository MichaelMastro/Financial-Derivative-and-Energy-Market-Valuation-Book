function [rootn, wn, Hcoef] =...
    HermiteWeightAndRoots (nend, OutputFlag  )
% HermiteWeightAndRoots calculates Coefficients of H_n
% The roots are calculated for each H_n polynomial. The Roots
% serve as abscissa points for Gauss-Hermite Quadrature
% Corresponding Gauss-Hermite Quadrature weigths at each root
% are also calculated. 

if (nargin <2), OutputFlag=1; end
if (nargin <1),  nend=12; end

set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold') 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);
close all

if (OutputFlag==1)
    figure; 
    nstart=1;
    HcoefPrev=HCoeffCalc(nstart-1);
    for n=nstart:nend
        Hcoef=HCoeffCalc(n);  
        %The roots of H_n polynomial returned in a column vector
        rootn=roots(Hcoef);
        wn=(2^(n-1)*factorial(n)*sqrt(pi) )./...
            (n.^2*(polyval(HcoefPrev,rootn).^2));
        HcoefPrev=Hcoef;
            fprintf(1, 'H(%g) Coefficients:',n);
            disp(Hcoef)
            fprintf(1, 'H(%g) Roots:',n);
            disp(rootn')    
            fprintf(1, 'H(%g) Weights:',n);
            disp(wn')
            fprintf(1, '\n');
        MFC=[n n/1.5 n/2]'/nend;
        plot(rootn,wn,'o','MarkerSize',...
            5*nend/n,'MarkerFaceColor',MFC)
        hold on  
    end
    hold off
    xlabel('Roots of H_n Hermite Polynomial'); 
    ylabel('Weights for H_n'); 
    title('Weights and Roots for Gauss-Hermite Quadrature')
    text(0.7,1.3,'Symbol Size ~ 1/n')
else % (OutputFlag==0) % No Output ; No Plot
    [Hcoef HcoefPrev] =HCoeffCalc(nend); 
    %The roots of H_n polynomial returned in a column vector
    rootn=roots(Hcoef);
    wn=(2^(nend-1)*factorial(nend)*sqrt(pi) )./...
        (nend.^2*(polyval(HcoefPrev,rootn).^2));
end

%weightsum=sum(wn) % Weights always add up to 1.7725
end

function [Hcoeff HcoeffPrev] = HCoeffCalc(n) 
% HCoeffCalc Calculates Hermite Polynomial Coefficients
% by starting with the known H(n=0) and H(=1) polynomial, 
% then iteratively calculating up to the H(n) polynomial
% H2(x)=4x^2-2      % p2 = [4 0 -2]
% H3(x)=8x^3-12x    % p3 = [8 0 -12 0]
Hcoeff =zeros(1,n+1);
if (n<1),Hcoeff(1)=1;  % n=0
elseif (n==1), Hcoeff(1)=2; % n=1
else
    Hcoeffnm1 = [1]; %Initialize H(n-1)=H(0)=1
    Hcoeffn= [2 0]; %Initialize H(n)= H(1) = 2x; 
    m=1;
    % H(n+1)=2xH(n)-2nH(n-1)
    while (m<n)       
        Hcoeffnp1=conv([2 0],Hcoeffn); % 2x*H(n)    
        SecondTerm=2*m*Hcoeffnm1; % 2n*H(n-1)
        len=length(SecondTerm); % find number of Coeffs to add    
        Hcoeffnp1(end-len+1:end)=Hcoeffnp1(end-len+1:end)-SecondTerm;
        Hcoeffnm1=Hcoeffn; % shift for next loop
        Hcoeffn=Hcoeffnp1; % shift for next loop
        m=m+1;
    end
    Hcoeff=Hcoeffnp1;
    HcoeffPrev=Hcoeffnm1;
end
end



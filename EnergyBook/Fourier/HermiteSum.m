function [] = HermiteSum(maxH)
% HermiteSum calculates exp(-x^2) by a Talyor summation
% with and without Hermite Polynomials 
% Both approaches are mathematically equivalent
% The third plot is first few Hermite polynomials over
% a range of x0 points. 

close all; clc
set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold') 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

if nargin ==0
    maxH=7;  
    % 4 terms acceptable over x=-0.1 to 0.1
    % 7 terms acceptable over x=-1 to 1 
    % Odd number of terms has better symmetry
    % since is actually even number of terms 
    % when H_0 is neglected
end
x=-1:.1:1;  
y=exp(-x.^2); % True value
x0=0.1;

% Calculate Taylor series without Hermite Polynomials
% Taylor = sum ( (x-a)^n/n! * dn/dxn(f)|@x0 )
% use first 4 terms
SimpleTay = zeros(maxH,length(x));
SimpleTay(1,:)= 1 ;
SimpleTay(2,:)= (x-x0).*     exp(-x0.^2).* (-2.*x0) / 1 ;
SimpleTay(3,:)= ((x-x0).^2).*exp(-x0.^2).*(4.*x0.^2-2) / 2;
SimpleTay(4,:)= ((x-x0).^3).*exp(-x0.^2).*(-8.*x0.^3+12.*x0)...
    / (3*2);
SimpleTaySum = sum(SimpleTay);
figure
plot(x,SimpleTaySum,':+',x,y,'--o',x,SimpleTay)
xlabel('x');  ylabel('y');
title('Taylor Series with 4 Terms')
legend('Taylor near x_0=0','y=e^{-x^2}', 'Location', 'South')

% Calculate Taylor series with maxH Hermite Polynomial terms
% Equivalent to previous calculation if maxH = 4
tay = zeros(maxH,length(x));
H = zeros(maxH,1);
% Recursive Calculation of Hermite Polynomials
H(1)=1;
H(2)=2*x0;
% Hn+1(x)=2xHn(x) - 2nHn-1(x)
for ind=2:(maxH-1)
    realn=ind-1;
    % real n is n-1; Matlab cannot index starting at 0
    H(ind+1)=2*x0*H(ind)-2*(realn).*H(ind-1);    
end

for ind=1:maxH
    tay(ind,:)=(((x-x0).^(ind-1)).* (-1)^(ind-1).*H(ind)...
                    .*exp(-x0.^2)) / factorial(ind-1);
end
taysum=sum(tay);

figure 
plot(x,taysum,':+',x,y,'--o',x,tay); 
legend('Taylor near x_0=0','y=e^{-x^2}', 'Location', 'South'); 
xlabel('x');  ylabel('y');
maxHstring = num2str(maxH);
title(['Summation via ' maxHstring ' Hermite Polynomials'])

% Use of Hermite polynomials above used only x0 point
% Plot first 4 Hermite polynomials over a range of x0 points
Herm=zeros(maxH,length(x));
Herm(1,:)=ones(1,length(x));
Herm(2,:)=2*x;
% Hn+1(x)=2xHn(x) - 2nHn-1(x)
for index=2:4
    realn=index-1;
    % real n is n-1; Matlab cannot index starting at 0
    Herm(index+1,:)=2*x.*Herm(index,:)-2*(realn).*Herm(index-1,:); 
end
figure
plot(x,Herm); legend('H_0', 'H_1', 'H_2', 'H_3', 'H_4')
xlabel('x_0'); ylabel('y'); title ('Hermite polynomials')
end


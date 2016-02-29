function [ Gab ] = GammaRand ( a,b )
% GammaRand generates Gamma Random variables 
% Pseudo-Code and description in Deville
% Levy Processes in Finance
% Requires a<=1
% usually simulate as Gamma (a1=a*deltaT,b)
% since deltaT usually small then a1 is small

% Since X/c ~= Gamma(a,b1=b*c) 
% calc. Gamma(a,1) 
% Gamma(a,b) ~= Gamma(a,1)  / b 

if (nargin <2), b=10; end
if (nargin <1), a=10; end

%if (a<=1) 
    % Ahrens-Dieter Gamma generator
    m=(a+exp(1))/exp(1);
    reject=1;
    while (reject==1)
        u1 = rand;    
        u2 = rand;
        Y=m*u1;
        if (Y <= 1)
            Z=Y^(1/a);
            if ( u2 < (exp (-Z)) )
                reject = 0; %accept = true
            end
        else
            Z=-log((m-Y)/a);
            if ( u2 <= (Z^(a-1)) )
                 reject = 0; %accept = true
            end
        end
    end    
Gab=Z/b; % = Gamma(a,1)/b
end


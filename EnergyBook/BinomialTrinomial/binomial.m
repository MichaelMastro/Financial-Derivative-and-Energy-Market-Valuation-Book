function OptionValue=binomial (S0,K,r,sigma,T,Steps);
%binomial calculates call price via recursive binomialgraph
%A futures option requires y=r thus b=0; 
% and option value equation without discounting
%An American option requires adding a premium 
% for early exercise

close all;
clc
set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold'); 
set(0,'defaultaxesfontsize',14); set(0,'defaulttextfontsize',14);

AmericanFlag=0; %1=American, other=European
FuturesFlag=0; %1=Futures
if nargin < 6
   S0=50;
   K=50;
   r=0.05;
   sigma=0.1;
   T=1;
   Steps=7;
    if (FuturesFlag==1)
        y=r; %Futures cost of carry b=0=r-y -> set y=r -> a=1
    else
        y=0; %yield % dividend or convenience yield; 
    end  
end


dt=T/Steps;
u=exp(sigma*sqrt(dt))
d=exp(-sigma*sqrt(dt))
a=exp((r-y)*dt)    %Growth Factor
p=(a-d)/(u-d)
%OneMinP=1-p
t0=0;
figure

OptionValue = binomialgraph(t0,S0,u,d,Steps,dt,K,p,r,...
                    AmericanFlag,FuturesFlag)
xlabel('Time (Years)'); ylabel('Asset and Call Value');
if (FuturesFlag==1)
    title(['Futures Call Binomial Tree: b = 0'])
else
    title(['Stock Call Binomial Tree: Yield = ',num2str(y)])
end

%
axis tight
hold off
end

function OptionValue = binomialgraph(t,S,u,d,Steps,dt,...
                            K,p,r,AmericanFlag,FuturesFlag)
% binomialgraph recursively steps through tree
% Start at expiration where OptionValue=max(0,(S-K))
% and recursively work back to present time t0 
% and price S0
% The Matlab Financial Toolbox has several built-in commands
% to calculate financial derivatives via a tree and to view 
% graphical tree. Since not a standard feature the following 
% code will graph the tree using basic Matlab plotting commands
% The recursive function is useful to generate graph and 
% recursion is intuitive to working option price backwards
% However, to simply calculate price without graph
% a 'lattice' tree converted to a vector should be faster 

   if (t<=((Steps-1)*dt)) 
       clear x y
       x(:,1)=[t (t+dt)];
       x(:,2)=x(:,1);
       y(:,1)=[S (S*u)];
       y(:,2)=[S (S*d)] ;   
% value of current node equals future node mult. by 
% prob. of up-movement p and future node time mult.
% by down-movement prob. 1-p   

         %European Value on asset discount at risk-free rate
     OptionValue=exp(-r*dt)*( p*binomialgraph(t+dt,S*u,u,d,...
                  Steps,dt,K,p,r,AmericanFlag,FuturesFlag)+...
      (1-p)*binomialgraph(t+dt,S*d,u,d,Steps,dt,K,p,r,...
                        AmericanFlag,FuturesFlag));

        %American call adds value of early exercise = (S-K)
         if (AmericanFlag==1)
                OptionValue=max((S-K),OptionValue); 
         end     
       plot(x,y)
       txstr1(1) = {num2str(S,'%.2f')};
       txstr1(2) = {num2str(OptionValue,'%.2f')}; 
       text(t,S,txstr1);
       hold on
   else  % Expiration nodes
       OptionValue=max(0,(S-K));
       txstr2(1) = {num2str(S,'%.2f')};
       txstr2(2) = {num2str(OptionValue,'%.2f')}; 
       text(t*0.95,S,txstr2);
       hold on
   end
end
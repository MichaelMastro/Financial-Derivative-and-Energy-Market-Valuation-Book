function Pz= signif (z);
lz=length(z);
j=1:1000; %assume 1000 ~ infinity for this calculation
for i=1:lz
    Pz(i)=2*sum((-1).^(j-1).*exp(-2.*j.^2.*z(i).^2));
end

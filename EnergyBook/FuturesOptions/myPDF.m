function pdf = myPDF (x, av, sg)
pdf = (1/sqrt(2*pi*sg*sg)).*exp( -((x-av).^2) / (2*sg*sg) );
end
%-------------------------------------------------------------
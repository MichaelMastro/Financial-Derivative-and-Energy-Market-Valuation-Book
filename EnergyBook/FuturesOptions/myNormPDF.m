function npdf = myNormPDF (x)
 
npdf = (1/sqrt(2*pi))*exp(-(x.*x)/2);
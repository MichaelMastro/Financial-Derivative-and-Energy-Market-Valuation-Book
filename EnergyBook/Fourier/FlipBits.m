%dr = bin2dec(fliplr(dec2bin(d,n))); 
N=16
d=0:N-1
n=log2(N)
bi=dec2bin(d)
bin=dec2bin(d,n)
flr=fliplr(dec2bin(d))
DecF=bin2dec(fliplr(dec2bin(d)))
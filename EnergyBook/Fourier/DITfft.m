function H = DITfft(data)
%DITfft graphs radix-2 decimation in time fft
%produces same numbers as Matlab fft

close all;
clc
set(0,'defaultaxeslinewidth',2); set(0,'defaultlinelinewidth',2);
set(0,'defaultaxesFontWeight','bold');set(0,'defaulttextFontWeight','bold') 
set(0,'defaultaxesfontsize',8); set(0,'defaulttextfontsize',14);

if (nargin == 0) 
    N=16;
   %use data of 0,1,2,...N-1
   %to make bin-reversal obvious on plot
   data=[0:N-1]; 
   %data=ones(1,N)%rand(1,N)
   %data(1)=0;
end

N=length (data(:))
MaxBin=log2(N)
%if data length not a power of two then truncate
%if want non power of 2 vectors then use Matlab's fft
if (MaxBin~=floor(MaxBin))
        MaxBin=floor(MaxBin);
        data=data(1:2^MaxBin);
        N=length (data(:))
end

%Matlab now has bitrevorder command in signal toolbox
%but interesting to see steps to flip binary index
DataIndex=0:N-1;
%bin=dec2bin(DataIndex,MaxBin);
%flr=fliplr(dec2bin(DataIndex,MaxBin));
FlipIndex=bin2dec(fliplr(dec2bin(DataIndex,MaxBin)));
H(FlipIndex+1)=data(DataIndex+1);
SaveFlip=H; % save for graphing

%Start at bottom (level =1) of tree with N vectors of length 1   
%and combine (y=even+W*odd; even-W*odd) until have
%one vector of length N
%
for level = 1:MaxBin
    StepSize=2^level
    HalfStep=StepSize/2;
    StepEnd=StepSize/2-1;
    
    k=[0:StepSize/2-1]
    w=exp(-i*2*pi.*k/StepSize) %array of twiddle factors   
     
    for Step = 1:StepSize:N-1 
    %Each step is an Even and Odd pair to be combined
%level:1 on bit reversed array
%   Step 1: Eindex=1, Oindex=2; Step 2: Eindex=3,Oindex=4;...
%level:2 take array for level 1
%   Step 1: Eindex=1 2, Oindex=3 4; Step 2: Eindex=5 6,Oindex=7 8;...
%level:2 take array for level 2
%   Step 1: Eindex=1 2 3 4, Oindex=5 6 7 8;
%continue until even and odd arrays length =N/2
     
        EvenIndex=(Step:(Step+StepEnd));  
        OddIndex=HalfStep+EvenIndex;
                   
        Feven=H(EvenIndex); Fodd=H(OddIndex);
%could combine 4 or 6 steps into 2 steps vs. emphasize symmetry
        H(EvenIndex)=Feven+w.*Fodd;
        H(OddIndex)=Feven-w.*Fodd;
    %End of fft math
        
    %Draw tree 
        Xmid=(EvenIndex(1)+OddIndex(1))/2;
        EvenX=[EvenIndex(1) Xmid]+level-2;
        OddX=[OddIndex(1) Xmid]+level-2;
        Y=[level^1.5 level^1.5+0.75];
        plot (EvenX, Y, OddX, Y); hold on
    %plot calculated data onto tree
        texstringMid(1)={(['E+W_n_kO'])};
        texstringMid(2)={num2str([H(EvenIndex)]','%5.1f')};
        texstringMid(3)={(['E-W_n_kO'])};
        texstringMid(4)={num2str([H(OddIndex)]','%5.1f')};
        text(EvenX(2)-level/3,Y(2)+0.5, texstringMid,...
            'EdgeColor','black','Margin',1)
    %%plot twiddle factors (W) onto odd branch
        OmegaStr(1)={(['W_n_k'])};
        OmegaStr(2)={num2str(w,'%5.1f')}; %
        text(EvenX(2)+level/5,Y(2)-0.25, OmegaStr,...
            'EdgeColor','red','LineStyle','--','Margin',1)
    end
end
for loop=0:N-1 % plot bin-reversed data at bottom of tree
    datastr=num2str(SaveFlip(loop+1));
    text(loop,0.75, datastr)   ;
end
xlim([-0.2 1.1*loop])
ylim([0.6 1.4*max(Y)])
set(gca, 'XTickLabelMode', 'Manual')
set(gca, 'XTick', [])

%figure %same as matlab fft within numerical error
%subplot(2,1,1);plot(real(H-fft(data)));%=0 
%subplot(2,1,2);plot(imag(H-fft(data)));%=0 
end

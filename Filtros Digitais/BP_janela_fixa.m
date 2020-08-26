%%Filtro digital FIR de janela fixa BP
clear all
close all
clc

%especificacoes
fa = 1000; %freq amostragem
f1 = 150; %fs1
f2 = 200; %fp1
f3 = 300; %fp2
f4 = 380; %fs2
Ap = 1; %atenuacao max na passagem
As = 40; %atenuacao min stop
Gtopo = 10; %ganho em 0 linear

fs1 = f1/(fa/2);
fp1 = f2/(fa/2);
fp2 = f3/(fa/2);
fs2 = f4/(fa/2);

wc1 = ((fs1+fp1)/2)+0.0065;
wc2 = ((fs2+fp2)/2)+0.005;

N = 51; %Ordem do filtro inicialmente estimada
if mod(N, 2) == 1  
    M = (N+1)/2; %numero de coeficientes se N impar
else 
    M = N/2; %numero de coeficientes de N par
end

CBP =((wc2*sinc(wc2*(-M:M)))-(wc1*sinc(wc1*(-M:M))));


RET = CBP;
HAMMING = CBP.*hamming(2*M+1)';
BLACKMAN = CBP.*blackman(2*M+1)';
BARTLETT = CBP.*bartlett(2*M+1)';
HANN = CBP.*hann(2*M+1)';

HAMMING = HAMMING*10^((Gtopo-0.5)/20);
HANN = HANN*10^((Gtopo-0.5)/20);
BARTLETT = BARTLETT*10^((Gtopo-0.5)/20);
BLACKMAN = BLACKMAN*10^((Gtopo-0.5)/20);

figure(1)
impz(HAMMING)
title('Resposta ao impulto do filtro BP de janela fixa de Hamming')

figure(2)
zplane(HAMMING)
title('Diagrama dos polos e zeros do filtro BP de Hamming')



[Hret, w] = freqz(HAMMING, 1);
figure(3)
subplot(211)
plot(w/pi, 20*log10(abs(Hret))); hold on;
title('Filtro BP de janela fixa HAMMING n=51')
 plot([fp1 fp1 fp2 fp2], [-(As+30) 9 9 -(As+30)], ':m')
 plot([0 fs1 fs1 fs2 fs2 1], [-30 -30 10 10 -30 -30], ':m')
 ylim([-(As+30) Ap/2+10])
 subplot(212)
 plot(w/pi, 20*log10(abs(Hret))); hold on;
 plot([fp1 fp1 fp2 fp2], [-(As+30) 9 9 -(As+30)], ':m')
 plot([0 fs1 fs1 fs2 fs2 1], [-30 -30 10 10 -30 -30], ':m')
 ylim([8 11])
 hold off;



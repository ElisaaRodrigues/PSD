%%Filtro digital FIR de janela fixa LP
clear all
close all
clc

%espeficicacoes
fa = 500;
f1 = 20; %freq passagem
f2 = 100; % freq stop
Ap = 3; %dB atenuacao max banda passagem
As = 35; %dB atenuacao min band rejeicao
Gtopo = 1; %dB ganho na banda de passagem linear

%dados calculados a partir
G0 = 0; %dB no prototipo ganho na banda de passagem
Gp = Gtopo-Ap;
Gs = Gtopo-As;

fp = f1/(fa/2);
fs = f2/(fa/2);
Wp = fp;
Ws = fs;
Wc = (Wp+Ws)/2;

N = 17; %Ordem do filtro inicialmente estimada
if mod(N, 2) == 1  
    M = (N+1)/2; %numero de coeficientes se N impar
else 
    M = N/2; %numero de coeficientes de N par
end

CLP = Wc*sinc(Wc*(-M:M));
RET = CLP;
HAMMING = CLP.*hamming(2*M+1)';
BLACKMAN = CLP.*blackman(2*M+1)';

HAMMING = HAMMING*10^((Gtopo-0.5)/20);
figure(1)
impz(HAMMING)
title('Resposta ao impulto do filtro LP de janela fixa')

figure(2)
zplane(HAMMING)
title('Diagrama dos polos e zeros do filtro LP Hamming')

[Hret, w] = freqz(HAMMING, 1);
figure(3)
plot(w/pi, 20*log10(abs(Hret))); hold on;
title('Filtro LP de janela fixa HAMMING n=17')
 plot([0 Wp Wp], [-2 -2 -(As+30)], ':m')
 plot([0 Ws Ws 1], [1 1 -34 -34], ':m')
 ylim([-(As+30) Ap/2+10])
 hold off;
 



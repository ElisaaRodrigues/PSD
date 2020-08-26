%%Filtro digital FIR de janela ajustavel BP
clear all
close all
clc

%especificacoes
fa = 1000; %freq amostragem
f1 = 150; % fs1
f2 = 200; %fp1
f3 = 300; %fp2
f4 = 380; %fs2

Ap = 1; %dB atenuacao max na passagem
As = 40; %dB atenucacao min na rejeicao
Gtopo = 10; %dB ganho em 0 linear

ws1 = f1/(fa/2);
ws2 = f4/(fa/2);
wp1 = f2/(fa/2);
wp2 = f3/(fa/2);


freq = [f1 f2 f3 f4]; %vetor de frequencias
a = [0 1 0]; %definindo filtro LP
dev = [10^(-As/20) 1-10^(-Ap/20) 10^(-As/20)]; %vetor de desvios (atenuacoes)

[n,wc,beta,ftype] = kaiserord(freq,a,dev,fa); %funcao para obter ordem e o beta 
n= n-6;
wc = [0.36 0.65]
wkaiser = kaiser(n+1,beta);

h_fir = fir1(n,wc,ftype,wkaiser,'noscale');

h_fir = h_fir*10^((Gtopo-0.3)/20);

figure(1)
impz(h_fir)
title('Resposta ao impulto do filtro BP de janela ajustavel de Kaiser')

figure(2)
zplane(h_fir)
title('Diagrama dos polos e zeros do filtro BP de Kaiser')

[Hw,w] =freqz(h_fir);
figure(3)
subplot(211)
plot(w/pi, 20*log10(abs(Hw))); hold on;
title(['Filtro BP de janela ajustavel KAISER n= ',num2str(n) ])
plot([wp1 wp1 wp2 wp2], [-(As+30) 9 9 -(As+30)], ':m')
plot([0 ws1 ws1 ws2 ws2 1], [-30 -30 10 10 -30 -30], ':m')
ylim([-(As+30) Ap/2+10])
subplot(212)
plot(w/pi, 20*log10(abs(Hw))); hold on;
plot([wp1 wp1 wp2 wp2], [-(As+30) 9 9 -(As+30)], ':m')
plot([0 ws1 ws1 ws2 ws2 1], [-30 -30 10 10 -30 -30], ':m')
ylim([8 11])




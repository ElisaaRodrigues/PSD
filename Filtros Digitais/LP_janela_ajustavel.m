%%Filtro digital FIR de janela ajustavel LP
clear all
close all
clc

%especificacoes
fa = 500; %freq amostragem
f1 = 20; % freq passagem
f2 = 100; %freq stop
Ap = 3; %dB atenuacao max na passagem
As = 35; %dB atenucacao min na rejeicao
Gtopo = 1; %dB ganho em 0 linear
ws = f1/(fa/2);
wp = f2/(fa/2);

freq = [f1 f2]; %vetor de frequencias
a = [1 0]; %definindo filtro LP
dev = [1-10^(-Ap/20) 10^(-As/20)]; %vetor de desvios (atenuacoes)

[n,wc,beta,ftype] = kaiserord(freq,a,dev,fa); %funcao para obter ordem e o beta 
n= n-4;
wc = wc-0.12;
wkaiser = kaiser(n+1,beta);

h_fir = fir1(n,wc,ftype,wkaiser,'noscale');

h_fir = h_fir*10^((Gtopo+2)/20);

figure(1)
impz(h_fir)
title('Resposta ao impulto do filtro LP de janela ajustavel')

figure(2)
zplane(h_fir)
title('Diagrama dos polos e zeros do filtro LP Kaiser')


[Hw,w] =freqz(h_fir);
figure(3)
plot(w/pi, 20*log10(abs(Hw))); hold on;
title(['Filtro LP de janela ajustavel KAISER n=',num2str(n) ])
plot([0 ws ws], [-2 -2 -(As+30)], ':m')
plot([0 wp wp 1], [1 1 -34 -34], ':m')
ylim([-(As+30) Ap/2+10])



%%Filtro HP com algoritmo de Park-McClellan
clear all
close all
clc

fa = 500; %freq amostragem
f1 = 20; %freq stop
f2 = 100; %freq passagem
Ap = 3; %ateanuacao max banda passagem
As = 35; %atenuacao min banda rejeicao
Gtopo = -20; %Ganho em 0 linear

freq = [f1 f2]; %definindo frequencias criticas
a = [0 1]; %definindo filtro HP
dev = [10^(-As/20) (10^(Ap/20)-1)/(10^(Ap/20)+1)]; %definindo ganhos max e min 'deltas'
[n,fo,ao,w_pm] = firpmord(freq,a,dev,fa);
n=10;
w_pm = [8 4]
b = firpm(n,fo,ao,w_pm);
bz = b*10^((Gtopo-0.5)/20);

figure(1)
impz(b)
title('Resposta ao impulto do filtro HP Parks-McClellan')

figure(2)
zplane(b)
title('Diagrama dos polos e zeros do filtro HP Parks-McClellan')

[h,w] = freqz(bz,1,1024,fa);
figure(3)
plot(w, 20*log10(abs(h))); hold on;
title('Filtro HP com algoritmo de Park-McClellan')
plot([f2 f2 fa/2], [-(As+30) -23 -23], ':m')
plot([0 f1 f1 fa/2], [-55 -55 Gtopo Gtopo], ':m')
ylim([-(As+30) Ap/2+10])
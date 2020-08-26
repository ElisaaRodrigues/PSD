%%Filtro LP com algoritmo de Park-McClellan
clear all
close all
clc

fa = 500; %freq amostragem
f1 = 20; %freq passagem
f2 = 100; %freq stop
Ap = 3; %ateanuacao max banda passagem
As = 35; %atenuacao min banda rejeicao
Gtopo = 1; %Ganho em 0 linear
Apa = Ap/2;

freq = [f1 f2]; %definindo frequencias criticas
a = [1 0]; %definindo filtro LP
dev = [(10^(Apa/20)-1)/(10^(Apa/20)+1)  10^(-As/20)]; %definindo ganhos max e min 'deltas'


[n,fo,ao,w] = firpmord(freq,a,dev,fa);
n=8;
b = firpm(n,fo,ao,w);

figure(1)
impz(b)
title('Resposta ao impulto do filtro LP Parks-McClellan')

figure(2)
zplane(b)
title('Diagrama dos polos e zeros do filtro LP Parks-McClellan')


[h,w] = freqz(b,1,1024,fa);
figure(3)
plot(w, 20*log10(abs(h))); hold on;
title('Filtro LP com algoritmo de Park-McClellan')
plot([0 f1 f1], [-2 -2 -(As+30)], ':m')
plot([0 f2 f2 fa/2], [1 1 -34 -34], ':m')
ylim([-(As+30) Ap/2+10])

%Alteracoes feitas na ordem do filtro e tambem na atenuacao da passagem na
%hora de projeta-lo
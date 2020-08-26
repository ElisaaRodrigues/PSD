%%Filtro BP com algoritmo de Park McClellan
clear all
close all
clc

fa = 1000;
f1 = 150;
f2 = 200;
f3 = 300;
f4 = 380;
Ap = 0.5;
As = 40;
Gtopo = 10;


freq = [f1 f2 f3 f4]; %definindo frequencias criticas
a = [0 1 0]; %definindo filtro HP
dev = [10^(-As/20) (10^(Ap/20)-1)/(10^(Ap/20)+1) 10^(-As/20)]; %definindo ganhos max e min 'deltas'
[n,fo,ao,w_pm] = firpmord(freq,a,dev,fa);
w_pm = [17 7 17]
n = 34;
b = firpm(n,fo,ao,w_pm);

bz = b*10^((Gtopo-0.5)/20);

figure(1)
impz(bz)
title('Resposta ao impulto do filtro BP Parks-McClellan')

figure(2)
zplane(bz)
title('Diagrama dos polos e zeros do filtro BP Parks-McClellan')




[h,w] = freqz(bz,1,1024,fa);
figure(3)
subplot(211)
plot(w, 20*log10(abs(h))); hold on;
title('Filtro BP com algoritmo de Park-McClellan')
plot([f2 f2 f3 f3], [-(As+30) 9 9 -(As+30)], ':m')
plot([0 f1 f1 f4 f4 fa/2], [-30 -30 10 10 -30 -30], ':m')
ylim([-(As+30) Ap/2+10])
subplot(212)
plot(w, 20*log10(abs(h))); hold on;
plot([f2 f2 f3 f3], [-(As+30) 9 9 -(As+30)], ':m')
plot([0 f1 f1 f4 f4 fa/2], [-30 -30 10 10 -30 -30], ':m')
ylim([8 11])





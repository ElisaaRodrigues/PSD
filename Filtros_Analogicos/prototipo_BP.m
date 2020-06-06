
%% Prototipo filtro DIGITAL BP  de Chebyshev Type II
clear all
close all
clc

%Especificacoes do filtro

fa = 1000; %frequencia de amostragem
f1 = 150; %fs1 frequencia de rejeicao 1
f2 = 200; %fp1 frequencia de passagem 1
f3 = 300; %fp2 frequencia de passagem 2
f4 = 380; %fs2 frequencia de rejeicao 2
Ap = 1; % atenuacao maxima na banda de passagem
As = 40; %atenuacao minima na banda de rejeicao
Gtopo = 10; %ganho em 0 na escala linear
G0 = 0; %ganho em 0 no prototipo
Gp = Gtopo-Ap; %Ganho na passagem em escala linear
Gs = Gtopo-As; %Ganho na rejeicao em escala linear

% Dados do filtro digital H(z)

teta_s1 = f1/(fa/2);
teta_s2 = f4/(fa/2);
teta_p1 = f2/(fa/2);
teta_p2 = f3/(fa/2);

% Dados do filtro analogico H(s)

lambda_s1 = 2*tan((teta_s1*pi)/2);
lambda_s2 = 2*tan((teta_s2*pi)/2);
lambda_p1 = 2*tan((teta_p1*pi)/2);
lambda_p2 = 2*tan((teta_p2*pi)/2);


% Especificacoes do Hp(p) LP prototipo

B = lambda_p2 - lambda_p1;
lambda0 = sqrt(lambda_p2*lambda_p1);

Omega_p = 1;

%Calculos abaixo para saber qual omega de rejeicao seria o mais restritivo
%Omega_s1 = abs((-lambda_s1^2+lambda0^2)/(B*lambda_s1))
%Omega_s2 = abs((-lambda_s2^2+lambda0^2)/(B*lambda_s2))
%tendo como resultado Omega_s1=2.2361 e Omega_s2=3.2774

Omega_s = abs((-lambda_s1^2+lambda0^2)/(B*lambda_s1));

% Determinando Hp(p) prototipo

[n,Omega_s] = cheb2ord(Omega_p,Omega_s,Ap,As,'s');
[b, a] = cheby2(n, As, Omega_s, 's');

[h, w] = freqs(b,a,10000);
mag = 20*log10(abs(h));
subplot(411)
plot(w, mag);
title('Hp(p) - LP prototipo')
ylim([-45 10]);
xlim([0 10]);
x1 = [0 1 1];
y1 = [-1 -1 -45];
line(x1, y1, 'Color', 'y', 'LineStyle', '-');
x2 = [0 1.8028 1.8028 10];
y2 = [0 0 -40 -40];
line(x2, y2, 'Color', 'y', 'LineStyle', '-');
grid on;


syms p
Np(p) = poly2sym(b, p);
Dp(p) = poly2sym(a, p);
Hp(p) = Np/Dp;
%pretty(vpa(Hp(p), 2))

syms s
p = ((s^2+lambda0^2)/(B*s));
Hs(s) = collect(Hp(p));
%pretty(vpa(Hs(s),2));

[Ns, Ds] = numden(Hs(s));
bs = sym2poly(Ns);
as = sym2poly(Ds);
%fazendo a normalizacao para obter coeficientes de menores valores e reduzir probabilidade de erros
bsn = bs/as(1);
asn = as/as(1);


[h2, w2] = freqs(bsn, asn, 10000);
mag2 = 20*log10(abs(h2));
subplot(412)
plot(w2, mag2)
title('Hs(s) - BP analogico')
ylim([-45 10]);
xlim([0 10]);
x1 = [1.4531 1.4531 2.7528 2.7528];
y1 = [-45 -1 -1 -45];
line(x1, y1, 'Color', 'g', 'LineStyle', '-');
x2 = [0 1.0191 1.0191 5.0514 5.0514 10];
y2 = [-40 -40 0 0 -40 -40];
line(x2, y2, 'Color', 'g', 'LineStyle', '-');
grid on;

%transformando para digital
syms z
s = 2*((z-1)/(z+1));
Hz(z) = collect(Hs(s));
%pretty(vpa(Hz(z), 2))
[Nz, Dz] = numden(Hz(z));
bz = sym2poly(Nz)*10^(Gtopo/20);
az = sym2poly(Dz);
%fazendo a normalizacao para obter coeficientes de menores valores e reduzir probabilidade de erros
bzn = bz/az(1);
azn = az/az(1);

[h3, w3] = freqz(bzn, azn, 10000);
mag3 = 20*log10(abs(h3));
subplot(413)
plot(w3/pi, mag3)
title('Hz(z) - BP digital')
ylim([-45 20]);
xlim([0 1]);
x1 = [0.4 0.4 0.6 0.6];
y1 = [-45 9 9 -45];
line(x1, y1, 'Color', 'm', 'LineStyle', '-');
x2 = [0 0.3 0.3 0.76 0.76 1];
y2 = [-30 -30 10 10 -30 -30];
line(x2, y2, 'Color', 'm', 'LineStyle', '-');
grid on;


subplot(414)
plot(w3/pi*(fa/2), mag3)
title('Hz(z) - BP digital com grafico na frequencia')
ylim([-45 20]);
xlim([0 500]);
x1 = [200 200 300 300];
y1 = [-45 9 9 -45];
line(x1, y1, 'Color', 'c', 'LineStyle', '-');
x2 = [0 150 150 380 380 500];
y2 = [-30 -30 10 10 -30 -30];
line(x2, y2, 'Color', 'c', 'LineStyle', '-');
grid on;










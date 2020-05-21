%% Prototipo filtro LP de Butterworth
% valores de teste sao do projeto
clear all
close all
clc 
G0 = 1;   %Ganho em 0 linear
A0 = 0;   %Ganho em 0 no prototipo
Ap = 3;   %Atenuacao na passagem
As = 25;   %Atenuacao em stop
Gp = G0 - Ap;   %Ganho na passagem
Gs = G0 - As;   %Ganho no stop
Wp = 2*pi*20;   %Frequencia de passagem
Ws = 2*pi*100;   %Frequencia de corte
Omega_p = Wp/Wp;   %Omega na passagem (W normalizado)
Omega_s = Ws/Wp;   %Omega no stop (W normalizado)

%E = sqrt((10^(0.1*Ap))-1);   %fator epsilon
E = 1;   %Para casos em que a atenuacao em Ap=3dB
n = (log10(10^(0.1*As)-1))/(2*log10(Omega_s));   %calculando ordem do filtro
n = ceil(n);    %funcao ceil arredonda para o proximo inteiro

%calculos para obter os polos dos filtros
p1 = exp(1j*pi*((2*1)+n-1)/(2*n));
p2 = exp(1j*pi*((2*2)+n-1)/(2*n));

%outra forma de calcular os polos:
k = 1:n;
pk = exp(1j*pi*((2*k)+n-1)/(2*n));
%zplane(1, poly(pk))

Dp = real(poly(pk)); %denominador da funcao de transferencia Hp (filtro prototipo)
%imag(poly(pk)) parte imaginaria mto pequena entao pode ser descartada

syms p 

%D(p) = (p - p1)*(p - p2);
Hp(p) = 1/(Dp(1)*p^2 + Dp(2)*p^1 + Dp(3)*p^0) %funcao de transferencia do prototipo Hp = 1/Dp
pretty(vpa(Hp(p),3))

[h, w] = freqs(1, Dp);
%semilogx(w, 20*log10(abs(h))); grid on; %plotando filtro prototipo
[h, w] = freqs(1, Dp, [0, 1, 5, 10]);
20*log10(abs(h)); %resultados das atenuacoes no prototipo

syms s 
Hs(s) = Hp(s/Wp) %obtendo a funcao de tranferencia do filtro linear
pretty(vpa(collect(Hs(s)), 7))

Ns = (15791.37)*10.^(G0/20); %Numerador da funcao de transferencia linear multiplicado pelo ganho linear em G0=4dB
Ds = [1 177.7153 15791.37]; %Denomidanor da funcao de trasnferencia linear
[h, w] = freqs(Ns, Ds);
semilogx(w, 20*log10(abs(h))); grid on; %plot do filtro linear
[h, w] = freqs(Ns, Ds, [0, Wp, Ws, 10*Ws]);
20*log10(abs(h)) %calculos dos ganhos no filtro linear

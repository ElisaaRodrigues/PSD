%% Prototipo filtro HP chebyshev 1
% valores sao do projeto 
clear all
close all
clc
f1 = 20;    %frequencia stop
f2 = 100;   %frequencia passagem
Ap = 2;     %Atenuacao na passagem
As = 60;    %Atenuacao no stop
A0 = 0;     %Atenuacao em 0 (G0 em dB)
G0 = -20;   %Ganho em 0 linear
Gp = G0-Ap;   %Ganho na passagem
Gs = G0-As;   %Ganho no stop
Ws = 2*pi*f1;   %frequencia de stop
Wp = 2*pi*f2;   %frequencia de passagem
Omega_p = 1;    %frequencia de passagem no prototipo
Omega_s = Wp/Ws;    %frequencia de stop no prototipo

E = sqrt((10.^(0.1*Ap))-1); %calculando valor do epsilon (Ap =! 3dB)
n = ceil((acosh (sqrt((10^(0.1*As)-1)/E^2)))/(acosh(Omega_s))); %calculando a ordem do filtro
k = 1:n;       %K varia de 1 a n
Phi2 = (1/n)*asinh(1/E); %calculos para obter os polos
Thetak = ((2*k-1)*pi)/(2*n); %calculos para obter os polos

pk = -1*sinh(Phi2)*sin(Thetak)+1j*cosh(Phi2)*cos(Thetak); %calculando os polos do filtro
Dp = real(poly(pk)); %vetor com o valores do numerados da funcao de transferencia Hp --- apenas a parte real de pk
%imag(poly(pk)) valores imaginarios pequenos, entao sao descartados
numerador_p = (Dp(end)/sqrt(1+(E^2))); %calculando o numerador do Hp, fazendo a multiplicacao pois d0 par

syms p
Hp(p) = numerador_p/(Dp(1)*p^4 + Dp(2)*p^3 + Dp(3)*p^2 + Dp(4)*p^1 + Dp(5)); %funcao de transferencia do prototipo
pretty(vpa(Hp(p), 3)); %funcao de transferencia do prototipo
[h, w] = freqs(numerador_p, Dp);
% figure(1)
% subplot(211)
% semilogx(w, 20*log10(abs(h))); grid on;
% subplot(212)
% semilogx(w, angle(h)); grid on;

[h, w] = freqs(numerador_p, Dp, [0, 1, 5, 10]); %calculando as atenuacoes no prototipo
20*log10(abs(h)) %resultados das atenuacao no prototipo

syms s
Hs(s) = Hp(Wp/s); %calculando a funcao de transferencia do filtro linear
pretty(vpa(collect(Hs(s)),7)) 

%[Ns, Ds] = numden(Hs(s)) 
%Valores de Ns e Ds obtidos a partir do pretty Hs(s), pela funcao numden
%vem a multiplicacao com o s tambem
Ns = [5.888728e15*10.^(G0/20) 0 0 0 0]; %numerador do Hs, multiplicado pelo ganho linear do filtro (G0 = -18dB)
Ds = [7.413469e15 1.169905e19 1.78717e22 6.400778e24 5.615252e27]; %denominador do Hs
[hs, ws] = freqs(Ns, Ds);
semilogx(ws, 20*log10(abs(hs))); grid on;
[hs, ws] = freqs(Ns, Ds, [0, Ws, Wp, 10*Wp]);

20*log10(abs(hs)) %resultados dos ganhos no filtro linear



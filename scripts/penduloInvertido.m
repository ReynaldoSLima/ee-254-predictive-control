function [A,B,C,D] = penduloInvertido(m, M, L, g, T)
%% Fun��o que retorna o modelo no espa�o de estados,
%  j� com as matrizes discretizadas no tempo
% Os par�metros de entrada s�o:
% [A,B,C,D] = penduloInvertido(m, M, L, g, T)
% m ----------- massa na extremidade do p�ndulo
% M ----------- massa do carrinho
% L ----------- comprimento da barra do p�ndulo
% g ----------- acelera��o da gravidade
% T ----------- per�odo de amostragem

%% Defini��es das matrizes no tempo cont�nuo:

Ac = [0 1 0 0;...
    0 0 m*g/M 0;...
    0 0 0 1;...
    0 0 (m/M+1)*g/L 0];

Bc = [0;...
     1/M;...
      0;...
     1/(M*L)];
 
Cc = [1 0 0 0;...
      0 0 1 0];

%% C�lculo das matrizes discretizadas, sa�das da fun��o  
  
A = expm(Ac*T);

syms tal;
B = double(int(expm(Ac*tal)*Bc,0,T));

C = Cc;

D = 0;
end
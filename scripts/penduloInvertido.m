function [A,B,C,D] = penduloInvertido(m, M, L, g, T)
%% Função que retorna o modelo no espaço de estados,
%  já com as matrizes discretizadas no tempo
% Os parâmetros de entrada são:
% [A,B,C,D] = penduloInvertido(m, M, L, g, T)
% m ----------- massa na extremidade do pêndulo
% M ----------- massa do carrinho
% L ----------- comprimento da barra do pêndulo
% g ----------- aceleração da gravidade
% T ----------- período de amostragem

%% Definições das matrizes no tempo contínuo:

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

%% Cálculo das matrizes discretizadas, saídas da função  
  
A = expm(Ac*T);

syms tal;
B = double(int(expm(Ac*tal)*Bc,0,T));

C = Cc;

D = 0;
end
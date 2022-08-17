function [K,Pf] = reguladorDLQR(A,B,Q,R)
%% Função que retorna o regulador DLQR,
%  já com as matrizes discretizadas no tempo
% Os parâmetros de entrada são:
% [K,Pf] = solucaoDLQR(A,B,Q,R)
% Do modelo xp = Ax + B discretizado:
% A ----------- Matriz do modelo discretizado A
% B ----------- Matriz do modelo discretizado B
% Q ----------- Matriz de pesos do estado
% R ----------- Matriz de pesos para o controle

[K,Pf,~] = dlqr(A,B,Q,R);
Af = A - B*K;
fprintf('Autovalores de A-BK:\n');
fprintf('%f\n',eig(Af));

end
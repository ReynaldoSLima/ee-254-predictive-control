function [r_min, r_max, So, bo] = limOInf(A, B, C, K, Sx, Su, bx, bu, x0)
%% Função que retorna os limites da referência dado o x0
% Os parâmetros de entrada são:
% [r_min, r_max, So, bo] = limOInf(A, B, C, K, Sx, Su, bx, bu, x0)
% Do modelo xp = Ax + B discretizado:
% A  ----------- Matriz do modelo discretizado A
% B  ----------- Matriz do modelo discretizado B
% C  ----------- Matriz do modelo discretizado C
% K  ----------- Matriz do regulador DLQR
% Sx ----------- Matriz de limites do estado Sx x <= bx
% bx ----------- Matriz de limites do estado Sx x <= bx
% Su ----------- Matriz de limites do controle Su u <= bu
% bu ----------- Matriz de limites do controle Su u <= bu
% x0  ----------- Estado inicial

ex = 1e-3*ones(size(Sx,1),1);
eu = 1e-3*ones(size(Su,1),1);
n = size(A,1);
q = size(C,1);
p = 1;
aux = [A-eye(n) B; C zeros(q,1)]\[zeros(n,q);eye(q)];
Nx = aux(1:n,:);
Nu = aux(n+1:end,:);
Af = A - B*K;

Achif = [Af B*(K*Nx + Nu);...
        zeros(q,n) eye(q)];

Gamma = [eye(n) zeros(n,q);...
         -K     (K*Nx + Nu);...
       zeros(n)     Nx;...
       zeros(p,n)   Nu];
   
Spsi = blkdiag(Sx,Su,Sx,Su);
    
bpsi = [bx;...
        bu;...
        bx - ex;...
        bu - eu];

max_iter = 100; tol = 1e-6;
[So, bo] = determina_oinf(Achif, Gamma, Spsi, bpsi, max_iter, tol);

Ox = So(:, 1:n);
Or = So(:, n+1:end);
Sf = Ox;
%bf = bo - Or*r;

br = bo - Ox*x0;

[row_p, col_p] = find(Or > 0);
[row_n, col_n] = find(Or < 0);
index_p1 = row_p(find(col_p == 1));
r_max1 = min(br(index_p1)./Or(index_p1,1));
index_p2 = row_p(find(col_p == 2));
r_max2 = min(br(index_p2)./Or(index_p2,2));
index_n1 = row_n(find(col_n == 1));
r_min1 = max(br(index_n1)./Or(index_n1,1));
index_n2 = row_n(find(col_n == 2));
r_min2 = max(br(index_n2)./Or(index_n2,2));

r_min = [r_min1; r_min2];
r_max = [r_max1; r_max2];

end
clear; close all; clc;
m = 0.1; M = 1; L = 0.2; g = 9.81; T = 0.1;

[A,B,C,D] = penduloInvertido(m, M, L, g, T);

n = size(A,1);
p = 1;
q = size(C,1);
Q = [10 0 0 0;...
      0 1 0 0;...
      0 0 50 0;...
      0 0 0 1];
R = 0.1;  

[K,Pf] = reguladorDLQR(A,B,Q,R);

Sx = [0 0 1 0;...
       0 0 -1 0];
    
bx = [pi/10;...
     pi/20];
 
Su = [1;...
     -1];
bu = [2;...
      2];
  
% A referência é em d
% Lembre que: x = [d; dp; th; thp]
% Então: r = [df; thf]
x0 = [0;0;0;0];
r = [2; 0]; % hehehe, fora de onde ele alcança
N = 15;

[r_min, r_max, So, bo] = limOInf(A, B, C, K, Sx, Su, bx, bu, x0);
Ox = So(:, 1:n);
Or = So(:, n+1:end);
Sf = Ox;
bf = bo - Or*r;
Oinf = Polyhedron('H',[Sf bf]);

kf = testeDM(A,B,C,K,Q,R,Pf,T,N,x0,r,m,M,g,L,Sx,bx,Su,bu,Sf,bf,Oinf);
testeDLQR(A,B,C,K,T,kf,x0,r,m,M,g,L,Sx,bx,Su,bu,Oinf);

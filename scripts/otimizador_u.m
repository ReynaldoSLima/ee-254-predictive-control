function un = otimizador_u(A,B,Q,R,Pf,N,Sx,bx,Su,bu,Sf,bf,x0,xbar,ubar)

n = size(A,1);
p = 1;
An = A; Rn = R; Sun = Su;
for i = 2:N
    An = [An;A^i];
    Rn = blkdiag(Rn,R);
    Sun = blkdiag(Sun,Su);
end
bun = repmat(bu,N,1);

Qn = Q; Sxn = Sx;
for i = 2:N-1
    Qn = blkdiag(Qn,Q);
    Sxn = blkdiag(Sxn,Sx);
end
Qn = blkdiag(Qn,Pf);
Sxn=blkdiag(Sxn,Sf);
bxn = [repmat(bx,N-1,1);bf];

Baux = mat2cell(zeros(n*N,p*N),n*ones(N,1),p*ones(1,N));

for i = 1:N
    for j = 1:i
        Baux{i,j} = (A^(i-j))*B;
    end
end
Bn = cell2mat(Baux);

Aqp = [Sxn*Bn;Sun];
bqp = [bxn - Sxn*An*x0;bun];
Hqp = Bn'*Qn*Bn + Rn;
Hqp = (Hqp + Hqp')/2;
fqp = (Bn'*Qn*(An*x0 - repmat(xbar,N,1)) - Rn*repmat(ubar,N,1));

un = quadprog(Hqp,fqp,Aqp,bqp);

end
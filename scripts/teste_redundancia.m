function t = teste_redundancia(S,b,c,d)
Sa = [S;c']; ba = [b;d+1];
x = linprog(-c,Sa,ba);
t = (c'*x - d);
end
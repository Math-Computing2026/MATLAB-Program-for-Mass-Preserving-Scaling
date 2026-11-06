







syms g1 g2 g3 g4 g5
a = sym((1:5)/5); 
A = [a.'.^1, a.'.^2, a.'.^3, a.'.^4, a.'.^5];
B = simplify(inv(A));
g = [g1; g2; g3; g4; g5];
x = simplify(B * g)
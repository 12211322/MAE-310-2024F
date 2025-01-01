function [xi, eta, w] = tri_Gauss2D(N1, N2)

n_all = find_tri_n_int(N1, N2);
xi = zeros(n_all);
eta = xi;
w = xi;

[x1, w1] = Gauss(N1, 0, 1);
[x2, w2] = Gauss(N2, 0, 1);

for ii = 1 :N1
    for jj = 1 : N2
        xi((jj - 1)*N1 + ii) = x1(ii);
        eta((jj - 1)*N1 + ii) = x2(jj);
        w((jj - 1)*N1 + ii) = w1(ii) * w2(jj);
    end
    N1 = N1 - 1;
    if N1 == 0
        N1 = 1;
    end
end
function [xi, eta, w] = tri_Gauss2D(N1, N2)

n_all = find_tri_n_int(N1, N2);
xi = zeros(n_all,1);
eta = xi;
w = xi;

[x1, w1] = Gauss(N1, 0, 1);
[x2, w2] = Gauss(N2, 0, 1);

% for ii = 1 :N1
%     for jj = 1 : N2
%         xi((jj - 1)*N1 + ii) = x1(ii);
%         eta((jj - 1)*N1 + ii) = x2(jj);
%         w((jj - 1)*N1 + ii) = w1(ii) * w2(jj);
%     end
%     N1 = N1 - 1;
%     if N1 == 0
%         N1 = 1;
%     end
% end

if N1 >= N2
    n_index = 0;
    n_here = N1;
    for ii = 1 : N2
        for jj = 1 : n_here 
            n_index = n_index + 1;
            xi(n_index) = x1(jj);
            eta(n_index) = x2(ii);
            w(n_index) = w1(jj) * w2(ii);
        end
        n_diff = ceil((n_here - 1) / (N2 - ii));
        n_here = n_here - n_diff;
        if ii == N2 - 1
            n_here = 1;
        end
    end
else
    n_index = 0;
    n_here = N2;
    for ii = 1 : N1
        for jj = 1 : n_here 
            n_index = n_index + 1;
            xi(n_index) = x1(ii);
            eta(n_index) = x2(jj);
            w(n_index) = w1(ii) * w2(jj);
        end
        n_diff = ceil((n_here - 1) / (N1 - ii));
        n_here = n_here - n_diff;
        if ii == N1 - 1
            n_here = 1;
        end
    end
end

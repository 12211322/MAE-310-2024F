function [n_int] = find_tri_n_int(N1, N2)

if N1 >= N2
    n_int = N1 - N2;
    n_int = n_int + N2 * (N2 + 1) / 2;
else
    n_int = N2 -N1;
    n_int = n_int + N1 * (N1 + 1) / 2;
end
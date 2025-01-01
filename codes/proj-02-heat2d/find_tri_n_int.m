function [n_int] = find_tri_n_int(N1, N2)

% if N1 >= N2
%     n_int = N1 - N2;
%     n_int = n_int + N2 * (N2 + 1) / 2;
% else
%     n_int = N2 - N1;
%     n_int = n_int + N1 * (N1 + 1) / 2;
% end

if N1 >= N2
    Nb = N1; Ns = N2;
else
    Nb = N2; Ns = N1;
end

n_int = 1 + Nb;

n_up = Nb; %means the real-time up-number about the calculation of the difference
n_low = 1; %it is always 1
for ii = 1 : Ns - 2
    n_diff = ceil((n_up - n_low) / (Ns - ii));
    n_up = n_up - n_diff;
    n_int = n_int + n_up;
end

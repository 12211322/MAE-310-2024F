clear all
%physical property
v = 0.3;
E = 1e9;
G = E/2/(1+v);

%boundary coundary
fx = 0; fy = 0;



%exact solution


%quadrature rule
n_int_xi = 3;
n_int_eta = 3;
n_int = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

%mesh generation
quarter;
n_en = 4; %number of nodes in each element
n_el = size(msh.QUADS,1); %number of all elements
n_np = msh.nbNod; %number of all nodes
x_coor = msh.POS(:,1);
y_coor = msh.POS(:,2);

%IEN array
IEN = msh.QUADS(:,1:4);

%ID array
ID = -1 .* ones(n_np,2);
for ii = 1:size(msh.LINES,1)
    if msh.LINES(ii,3) == 11
        ID(msh.LINES(ii,1),2) = 0;
        ID(msh.LINES(ii,2),2) = 0;
    elseif msh.LINES(ii,3) == 9
        ID(msh.LINES(ii,1),1) = 0;
        ID(msh.LINES(ii,2),1) = 0;
    else
    end       
end
index = 0;
for ii = 1:n_np
    for jj = 1:2
        if ID(ii,jj) == -1
            index = index + 1;
            ID(ii,jj) = index;
        end
    end
end
n_eq = index;

%LM array
LM = zeros(n_el, 2*n_en);
for ii = 1:n_el
    for jj = 1:2*n_en
        LM(ii,jj) = ID(IEN(ii,ceil(jj/2)),rem(jj-1,2) + 1); %ceil为向上取整，rem为取余
    end
end

D = (E/(1-v*v)).*[1, v, 0;
    v, 1, 0;
    0, 0, 0.5-0.5*v]; %in 2D case

K = zeros(n_eq,n_eq);
F = zeros(n_eq,1);

%loop over element to assembly the matrix and vector
for ee = 1:n_el
    x_ele = x_coor(IEN(ee, :));
    y_ele = y_coor(IEN(ee, :));

    k_ele = zeros(2*n_en,2*n_en);
    f_ele = zeros(2*n_en,1);

    for ll = 1:n_int
        x_l = 0.0; y_l = 0.0;
        dx_dxi = 0.0; dx_deta = 0.0;
        dy_dxi = 0.0; dy_deta = 0.0;
        for aa = 1:n_en
            x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
            y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
            dx_deta = dx_deta + x_ele(aa) * Na_eta;
            dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
            dy_deta = dy_deta + y_ele(aa) * Na_eta;
        end

        detJ = dx_dxi * dy_deta - dx_deta * dy_dxi; %the Jacobian determinant

        for aa = 1 : n_en
            Na = Quad(aa, xi(ll), eta(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
            Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;

            f_ele(2*aa-1) = f_ele(2*aa-1) + weight(ll) * detJ * fx(x_l, y_l) * Na;
            f_ele(2*aa) = f_ele(2*aa) + weight(ll) * detJ * fy(x_l, y_l) * Na;

            for bb = 1 : n_en
                Nb = Quad(bb, xi(ll), eta(ll));
                [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
                Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
                Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;

                B_a = [Na_x, 0; 0, Na_y; Na_y, Na_x]; B_b = [Nb_x, 0; 0, Nb_y; Nb_y, Nb_x];
                k_ele(2*aa-1:2*aa,2*bb-1:2*bb) = k_ele(2*aa-1:2*aa,2*bb-1:2*bb) + (weight(ll) * detJ) .* B_a' * D * B_b;


            end % end of bb loop
        end % end of aa loop
    end % end of quadrature loop

    for ii = 1:2
        for aa = 1 : n_en
            pp = 2*(aa-1) + ii;
            PP = LM(ee, pp);
            if PP > 0
                F(PP) = F(PP) + f_ele(pp);
                for jj = 1:2
                    for bb = 1 : n_en
                        qq = 2*(bb-1) + jj;
                        QQ = LM(ee, qq);
                        if QQ > 0
                            K(PP, QQ) = K(PP, QQ) + k_ele(pp, qq);
                        else
                            % modify F with the boundary data
                            % here we do nothing because the boundary data g is zero or
                            % homogeneous
                        end
                    end
                end
            end
        end
    end
end

dn = K \ F;
x_disp = zeros(n_np, 1);
y_disp = zeros(n_np, 1);
for ii = 1 : n_np
    x_index = ID(ii,1);
    y_index = ID(ii,2);
    if x_index > 0
        x_disp(ii) = dn(x_index);
        if y_index > 0
            y_disp(ii) = dn(y_index);
        else
        end
    else
        % modify disp with the g data. Here it does nothing because g is zero
    end
end











function [output] = whole_rectangle(i,type,assum)
%i means the mesh sparseness; 
%type determinate the result, 1-error, 2-displacement, 3-strain, 4-stress
%assum means the different plane assumption, it will change the D

% physical property
E = 1e9; 
v = 0.3;
switch assum
    case 1
        D = (E/(1-v*v)).*[1, v, 0; v, 1, 0; 0, 0, 0.5-0.5*v]; %in 2D case
    case 2
        D = (E/(v+1)/(2*v-1)).*[1, v, 0; v, 1, 0; 0, 0, 0.5-0.5*v];
end

%exact solution
fx = @(x,y) (2*E*y*(y - 1))/(v^2 - 1) - (E*(v/2 - 1/2)*((x - 1)*(y - 1) + x*y + 2*x*(x - 1) + x*(y - 1) + y*(x - 1)))/(v^2 - 1) + (E*v*((x - 1)*(y - 1) + x*y + x*(y - 1) + y*(x - 1)))/(v^2 - 1);
fy = @(x,y) (2*E*x*(x - 1))/(v^2 - 1) - (E*(v/2 - 1/2)*((x - 1)*(y - 1) + x*y + x*(y - 1) + y*(x - 1) + 2*y*(y - 1)))/(v^2 - 1) + (E*v*((x - 1)*(y - 1) + x*y + x*(y - 1) + y*(x - 1)))/(v^2 - 1);
exactX = @(x,y) x*(1-x)*y*(1-y); exactY = @(x,y) x*(1-x)*y*(1-y);
exactX_x = @(x,y) (1-2*x)*y*(1-y); exactY_x = @(x,y) (1-2*x)*y*(1-y);
exactX_y = @(x,y) x*(1-x)*(1-2*y); exactY_y = @(x,y) x*(1-x)*(1-2*y); 
exactX_xx = @(x,y) -2*y*(1-y); exactY_xx = @(x,y) -2*y*(1-y);
exactX_xy = @(x,y) (1-2*x)*(1-2*y); exactY_xy = @(x,y) (1-2*x)*(1-2*y);
exactX_yy = @(x,y) -2*x*(1-x); exactY_yy = @(x,y) -2*x*(1-x);

% quadrature rule
n_int_xi  = 3; %take three quadrature points on ξ direction
n_int_eta = 3; %take three quadrature points on η direction
n_int     = n_int_xi * n_int_eta; %the number of all quadrature points in a 2D element
% here use arrays to remember 2D points' coordinates in a reference element and 2D points in this array's order is ξ1η1 ξ2η1 ... ξ1η2 ξ2η2 ...
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta); %where xi remember 2D points' ξcoordinates and weight's members are two 1D weights' product

% mesh generation
i = i * 20;
n_en   = 4;               % number of nodes in an element
n_el_x = i;               % number of elements in x-dir
n_el_y = i;               % number of elements in y-dir
qua_n_el   = n_el_x * n_el_y; % total number of elements

n_np_x = n_el_x + 1;      % number of nodal points in x-dir
n_np_y = n_el_y + 1;      % number of nodal points in y-dir
n_np   = n_np_x * n_np_y; % total number of nodal points

x_coor = zeros(n_np, 1);
y_coor = x_coor;

hx = 1.0 / n_el_x;        % mesh size in x-dir
hy = 1.0 / n_el_y;        % mesh size in y-dir

% generate the nodal coordinates
for ny = 1 : n_np_y
  for nx = 1 : n_np_x
    index = (ny-1)*n_np_x + nx; % nodal index
    x_coor(index) = (nx-1) * hx;
    y_coor(index) = (ny-1) * hy;
  end
end

% IEN array
qua_IEN = zeros(qua_n_el, n_en);
for ex = 1 : n_el_x
  for ey = 1 : n_el_y
    ee = (ey-1) * n_el_x + ex; % element index
    qua_IEN(ee, 1) = (ey-1) * n_np_x + ex;
    qua_IEN(ee, 2) = (ey-1) * n_np_x + ex + 1;
    qua_IEN(ee, 3) =  ey    * n_np_x + ex + 1;
    qua_IEN(ee, 4) =  ey    * n_np_x + ex;
  end
end

% ID array 
ID = zeros(n_np,2);
counter = 0;
for ny = 2 : n_np_y - 1 %ignore the first row and the last row
  for nx = 2 : n_np_x - 1 %ignore the first column and the last column
    index = (ny-1)*n_np_x + nx;
    counter = counter + 1;
    ID(index,1) = counter;
    counter = counter + 1;
    ID(index,2) = counter;
  end
end

n_eq = counter; %the number of equations in Kb = f, i.e. the number of unknow points

LM = zeros(qua_n_el, 2*n_en);
for ii = 1:qua_n_el
    for jj = 1:2*n_en
        LM(ii,jj) = ID(qua_IEN(ii,ceil(jj/2)),rem(jj-1,2) + 1); %ceil为向上取整，rem为取余
    end
end

K = zeros(n_eq,n_eq);
F = zeros(n_eq,1);

%loop over element to assembly the matrix and vector
for ee = 1:qua_n_el
    x_ele = x_coor(qua_IEN(ee, :));
    y_ele = y_coor(qua_IEN(ee, :));

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

%solve the displacement vector
dn = K \ F;

x_disp = zeros(n_np, 1);
y_disp = zeros(n_np, 1);
for ii = 1 : n_np
    x_index = ID(ii,1);
    y_index = ID(ii,2);
    if x_index > 0
        x_disp(ii) = dn(x_index);
    else
    end
    if y_index > 0
        y_disp(ii) = dn(y_index);
    else
    end
    
        % modify disp with the g data. Here it does nothing because g is zero
    
end

%output part
switch type
    case 1
        %error part
        e0 = 0.0; e1 = 0.0; u2 = 0.0;
        for ee = 1 : qua_n_el
            x_ele = x_coor(qua_IEN(ee, :));
            y_ele = y_coor(qua_IEN(ee, :));
            cx_ele = x_disp(qua_IEN(ee, :));
            cy_ele = y_disp(qua_IEN(ee, :));

            for ll = 1 : n_int
                x_l = 0.0; y_l = 0.0;
                dx_dxi = 0.0; dx_deta = 0.0; dy_dxi = 0.0; dy_deta = 0.0;
                uh_l = 0.0; uh_x_l = 0.0; uh_y_l = 0.0;

                for aa = 1 : n_en %this loop is to calculate the value of dx_dxi...
                    x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll)); %since the property of the node exactly, this x_l means the ture x coordinate of the integral point
                    y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));
                    [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
                    dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
                    dx_deta = dx_deta + x_ele(aa) * Na_eta;
                    dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
                    dy_deta = dy_deta + y_ele(aa) * Na_eta;
                end
                detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;

                %gain the coordinates of integral points, we can calculate uh ...
                u_l = exactX(x_l, y_l);
                u_x_l = exactX_x(x_l, y_l); u_y_l = exactX_y(x_l, y_l);
                u_xx_l = exactX_xx(x_l, y_l); u_xy_l = exactX_xy(x_l, y_l); u_yy_l = exactX_yy(x_l, y_l);

                for aa = 1 : n_en %this loop is to calculate the value of uh_l ...
                    c_ele = cx_ele(aa);
                    uh_l = uh_l + c_ele * Quad(aa, xi(ll), eta(ll));
                    [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));

                    uh_x_l = uh_x_l + c_ele * (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
                    uh_y_l = uh_y_l + c_ele * (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
                end

                %we should add weights in the integral loop
                e0 = e0 + weight(ll) * ((u_l - uh_l).^2);
                e1 = e1 + weight(ll) * ((u_l - uh_l).^2 + (u_x_l - uh_x_l).^2 + (u_y_l - uh_y_l).^2);
                u2 = u2 + weight(ll) * (u_l.^2 + u_x_l.^2 + u_y_l.^2 + u_xx_l.^2 + 2*(u_xy_l.^2) + u_yy_l.^2);
            end
        end

        L2_space = sqrt(e0./u2);
        H1_space = sqrt(e1./u2);
        output = [log(H1_space), log(L2_space), log(hx)];
    case 2
        %draw displacement
        figure;
        hold on;
        [X, Y] = meshgrid(0 : hx : 1, 0 : hy : 1);
        Z = reshape(x_disp, n_np_x, n_np_y)';
        surf(X, Y, Z);
        title('approximate x-displacement');
        xlabel('x');
        ylabel('y');
        shading interp;
        axis equal;
        colormap jet;
        colorbar;

        figure;
        hold on;
        [X, Y] = meshgrid(0 : hx : 1, 0 : hy : 1);
        Z = reshape(y_disp, n_np_x, n_np_y)';
        surf(X, Y, Z);
        title('approximate y-displacement');
        xlabel('x');
        ylabel('y');
        shading interp;
        axis equal;
        colormap jet;
        colorbar;
    case 3
end
end


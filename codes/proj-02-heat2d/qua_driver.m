kappa = 1.0; % conductivity

% exact solution
exact = @(x,y) x*(1-x)*y*(1-y);
exact_x = @(x,y) (1-2*x)*y*(1-y);
exact_y = @(x,y) x*(1-x)*(1-2*y);
exact_xx = @(x,y) -2*y*(1-y);
exact_xy = @(x,y) (1-2*x)*(1-2*y);
exact_yy = @(x,y) -2*x*(1-x);

% take_x = 0:0.01:1;
% take_y = 0:0.01:1;
% [xx, yy] = meshgrid(take_x, take_y);
% zz = xx.*(1-xx).*yy.*(1-yy);
% mesh(xx,yy,zz);

f = @(x,y) 2.0*kappa*x*(1-x) + 2.0*kappa*y*(1-y); % source term

all_data = zeros(5,3);

for i = 20:20:100

% quadrature rule
n_int_xi  = 3; %take three quadrature points on ξ direction
n_int_eta = 3; %take three quadrature points on η direction
n_int     = n_int_xi * n_int_eta; %the number of all quadrature points in a 2D element
% here use arrays to remember 2D points' coordinates in a reference element and 2D points in this array's order is ξ1η1 ξ2η1 ... ξ1η2 ξ2η2 ...
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta); %where xi remember 2D points' ξcoordinates and weight's members are two 1D weights' product

% mesh generation
n_en   = 4 ;               % number of nodes in an element
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
ID = zeros(n_np,1);
counter = 0;
for ny = 2 : n_np_y - 1 %ignore the first row
  for nx = 2 : n_np_x - 1 %ignore the first column
    index = (ny-1)*n_np_x + nx;
    counter = counter + 1;
    ID(index) = counter;  
  end
end

n_eq = counter; %the number of equations in Kb = f, i.e. the number of unknow points

LM = ID(qua_IEN);

% allocate the stiffness matrix and load vector
K = spalloc(n_eq, n_eq, 9 * n_eq);
F = zeros(n_eq, 1);

% loop over element to assembly the matrix and vector
for ee = 1 : qua_n_el
  x_ele = x_coor( qua_IEN(ee, 1:n_en) );
  y_ele = y_coor( qua_IEN(ee, 1:n_en) );
  
  k_ele = zeros(n_en, n_en); % element stiffness matrix
  f_ele = zeros(n_en, 1);    % element load vector
  
  for ll = 1 : n_int
    x_l = 0.0; y_l = 0.0;
    dx_dxi = 0.0; dx_deta = 0.0;
    dy_dxi = 0.0; dy_deta = 0.0;
    for aa = 1 : n_en %get each quadrature point's x = xieNi; take (-1,-1) (1,-1) (1,1) (-1,1) when aa = 1,2,3,4
      %Quad means the 1-order shape function on triangle
      x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll)); %since the property of the node exactly, this x_l means the ture x coordinate
      y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));    
      %Quad_grad means the 1-order shape function's derivative on triangle
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
      
      f_ele(aa) = f_ele(aa) + weight(ll) * detJ * f(x_l, y_l) * Na;
      
      for bb = 1 : n_en
        Nb = Quad(bb, xi(ll), eta(ll));
        [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
        Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
        Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;
        
        k_ele(aa, bb) = k_ele(aa,bb) + weight(ll) * detJ * kappa * (Na_x * Nb_x + Na_y * Nb_y);
      end % end of bb loop
    end % end of aa loop
  end % end of quadrature loop
 
  for aa = 1 : n_en
    PP = LM(ee, aa);
    if PP > 0
      F(PP) = F(PP) + f_ele(aa);
      
      for bb = 1 : n_en
        QQ = LM(ee, bb);
        if QQ > 0
          K(PP, QQ) = K(PP, QQ) + k_ele(aa, bb);
        else
          % modify F with the boundary data
          % here we do nothing because the boundary data g is zero or
          % homogeneous
        end
      end  
    end
  end
end

% solve the stiffness matrix
qua_dn = K \ F;

% insert dn back into the vector for all nodes
qua_disp = zeros(n_np, 1);

for ii = 1 : n_np
  index = ID(ii);
  if index > 0
    qua_disp(ii) = qua_dn(index);
  else
    % modify disp with the g data. Here it does nothing because g is zero
  end
end

%integral part
e0 = 0.0; e1 = 0.0; u2 = 0.0;
for ee = 1 : qua_n_el
    x_ele = x_coor(qua_IEN(ee, :));
    y_ele = y_coor(qua_IEN(ee, :));
    d_ele = qua_disp(qua_IEN(ee, :));

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
        u_l = exact(x_l, y_l);
        u_x_l = exact_x(x_l, y_l); u_y_l = exact_y(x_l, y_l);
        u_xx_l = exact_xx(x_l, y_l); u_xy_l = exact_xy(x_l, y_l); u_yy_l = exact_yy(x_l, y_l);

        for aa = 1 : n_en %this loop is to calculate the value of uh_l ...
            uh_l = uh_l + d_ele(aa) * Quad(aa, xi(ll), eta(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));

            uh_x_l = uh_x_l + d_ele(aa) * (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
            uh_y_l = uh_y_l + d_ele(aa) * (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
        end

        %we should add weights in the integral loop
        e0 = e0 + weight(ll) * ((u_l - uh_l).^2);
        e1 = e1 + weight(ll) * ((u_l - uh_l).^2 + (u_x_l - uh_x_l).^2 + (u_y_l - uh_y_l).^2);
        u2 = u2 + weight(ll) * (u_l.^2 + u_x_l.^2 + u_y_l.^2 + u_xx_l.^2 + 2*(u_xy_l.^2) + u_yy_l.^2);
    end
end

L2_space = sqrt(e0./u2);
H1_space = sqrt(e1./u2);


all_data(i/20,:) = [log(hx), log(L2_space), log(H1_space)];
end

plot(all_data(:,1),all_data(:,2));
hold on
plot(all_data(:,1),all_data(:,3));
legend('e in L2','e in H1');
% save the solution vector and number of elements to disp with name
% HEAT.mat
% save("HEAT", "disp", "n_el_x", "n_el_y");

% EOF
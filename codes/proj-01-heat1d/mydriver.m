% clear all; clc; clf; % clean the memory, screen, and figure

% Problem definition
f = @(x) -20*x.^3; % f(x) is the source
g = 1.0;           % u    = g  at x = 1
h = 0.0;           % -u,x = h  at x = 0

all_er = zeros(7,3);

for i = 2:2:16 %differnet number of physical elements, i.e. n_el
% for j = 1:6
% Setup the mesh
pp   = 3;              % shape function's polynomial degree
n_en = pp + 1;         % number of nodes in each element i.e. number of shape functions
n_el = i;              % number of physical elements
n_np = n_el * pp + 1;  % number of nodes in all physical elements
n_eq = n_np - 1;       % number of equations i.e. row# of F-vector
n_int = 10;            % number of integral sampling points in each element

hh = 1.0 / (n_np - 1); % space between two adjacent nodes
x_coor = 0 : hh : 1;   % nodal coordinates in all physical elements for equally spaced nodes

%Set up and finish IEN array
IEN = zeros(n_el, n_en);
for ee = 1 : n_el
  for aa = 1 : n_en
    IEN(ee, aa) = (ee - 1) * pp + aa;
  end
end

% Setup the ID array for the problem
ID = 1 : n_np;
ID(end) = 0;

% Set up the sampling points' coordinates and weights in reference element
[xi, weight] = Gauss(n_int, -1, 1);

% allocate the stiffness matrix
K = spalloc(n_eq, n_eq, (2*pp+1)*n_eq); %spalloc matrix can save memory space
F = zeros(n_eq, 1);

% Set up the stiffness matrix and load vector
for ee = 1 : n_el
  k_ele = zeros(n_en, n_en); % allocate a zero element stiffness matrix
  f_ele = zeros(n_en, 1);    % allocate a zero element load vector

  x_ele = x_coor(IEN(ee,:)); % take out each physical element's node coordinates

  % quadrature loop
  for qua = 1 : n_int %go through all intergal sampling points in element
    dx_dxi = 0.0;
    x_l = 0.0;
    for aa = 1 : n_en
      x_l    = x_l    + x_ele(aa) * PolyShape(pp, aa, xi(qua), 0); %get I.S.P coordinates of P.E. by mapping from R.E.
      dx_dxi = dx_dxi + x_ele(aa) * PolyShape(pp, aa, xi(qua), 1); %get dx/dξ at this I.S.P.
    end
    dxi_dx = 1.0 / dx_dxi; %get dξ/dx

    for aa = 1 : n_en %go through all shape funtions in R.E.
      f_ele(aa) = f_ele(aa) + weight(qua) * PolyShape(pp, aa, xi(qua), 0) * f(x_l) * dx_dxi;
      for bb = 1 : n_en %go through all shape funtions in R.E.
        k_ele(aa, bb) = k_ele(aa, bb) + weight(qua) * PolyShape(pp, aa, xi(qua), 1) * PolyShape(pp, bb, xi(qua), 1) * dxi_dx;
      end
    end
  end
 
  % Assembly of the matrix and vector based on the ID or LM data
  for aa = 1 : n_en
    P = ID(IEN(ee,aa));
    if(P > 0)
      F(P) = F(P) + f_ele(aa);
      for bb = 1 : n_en
        Q = ID(IEN(ee,bb));
        if(Q > 0)
          K(P, Q) = K(P, Q) + k_ele(aa, bb);
        else
          F(P) = F(P) - k_ele(aa, bb) * g; % handles the Dirichlet boundary data
        end
      end
    end
  end
end

% ee = 1 F = NA(0)xh
F(ID(IEN(1,1))) = F(ID(IEN(1,1))) + h;

% Solve Kd = F equation
d_temp = K \ F;

disp = [d_temp; g];

% Postprocessing: visualization
%plot(x_coor, disp, '--r','LineWidth',3);

%x_sam = 0 : 0.01 : 1;
%y_sam = x_sam.^5;
%hold on;
%plot(x_sam, y_sam, '-k', 'LineWidth', 3);

n_sam = 20;
xi_sam = -1 : (2/n_sam) : 1; %take 20 draw sampling points' coordinates in reference element

x_sam = zeros(n_el * n_sam + 1, 1);
y_sam = x_sam; % store the exact solution value at sampling points
u_sam = x_sam; % store the numerical solution value at sampling pts
ux_sam = x_sam;
yx_sam = x_sam;

erLU = 0;
erLL = 0;
erHU = 0;
erHL = 0;

for ee = 1 : n_el %go through all physical elements
  x_ele = x_coor( IEN(ee, :) ); %get eeth P.E.'s nodes' coordinates
  u_ele = disp( IEN(ee, :) ); %get eeth P.E.'s nodes' coefficients

  if ee == n_el
    n_sam_end = n_sam+1;
  else
    n_sam_end = n_sam;
  end

  for ll = 1 : n_sam_end %go through all intergal samping points in each P.E.
    x_l = 0.0;
    u_l = 0.0;
    ux_l = 0.0;
    for aa = 1 : n_en %go through all nodes in each P.E.
      x_l = x_l + x_ele(aa) * PolyShape(pp, aa, xi_sam(ll), 0); %get plot's x-axis point
      u_l = u_l + u_ele(aa) * PolyShape(pp, aa, xi_sam(ll), 0); %get uh_plot's y-axis point
      ux_l = ux_l + 2 * n_el * u_ele(aa) * PolyShape(pp, aa, xi_sam(ll), 1); %get uhx_plot's y-axis point
      
    end

    %store all needed statistic
    x_sam( (ee-1)*n_sam + ll ) = x_l;
    u_sam( (ee-1)*n_sam + ll ) = u_l;
    ux_sam( (ee-1)*n_sam + ll ) = ux_l;
    y_sam( (ee-1)*n_sam + ll ) = x_l^5;
    yx_sam( (ee-1)*n_sam + ll ) = 5*(x_l^4);

    erLL = erLL + (x_l^5).^2;
    erLU = erLU + (u_l - x_l^5).^2;
    erHL = erHL + (5*(x_l^4)).^2;
    erHU = erHU + (ux_l - 5*(x_l^4)).^2;

  end
end

erL = sqrt(erLU./erLL);
erH = sqrt(erHU./erHL);
all_er(i/2,:) = [log(erL),log(erH),log(hh)];

% plot(x_sam,ux_sam); %plot the slope
% hold on
% plot(x_sam,yx_sam);

% plot(x_sam, u_sam, '-r','LineWidth',3); %plot the function
% hold on;
% plot(x_sam, y_sam, '-k','LineWidth',1);
end

x = all_er(:,3);
y1 = all_er(:,1);
y2 = all_er(:,2);
plot(all_er(:,3),all_er(:,2));





























% EOF
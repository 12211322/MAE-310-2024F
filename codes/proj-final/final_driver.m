clear all
v = 0.3;
E = 1e9;
G = E/2/(1+v);

%exact solution

%define element shape
e_shape = 4;

%quadrature rule
n_int_xi = 3;
n_int_eta = 3;
if e_shape == 3
    n_int = find_tri_n_int(n_int_xi, n_int_eta);
    [xi, eta, weight] = tri_Gauss2D(n_int_xi, n_int_eta);
elseif e_shape == 4
    n_int = n_int_xi * n_int_eta;
    [xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);
else
    error('Error: e_shape must be 3 or 4.');
end

%mesh generation
n_en = e_shape;
if e_shape == 3
    tri_mesh;
    n_el = size(msh.TRIANGLES,1);
    col_index = [1, 2, 3];
    IEN = msh.TRIANGLES(:,col_index);
elseif e_shape == 4
    quad_mesh;
    n_el = size(msh.QUADS,1);
    col_index = [1, 2, 3, 4];
    IEN = msh.QUADS(:,col_index);
else
    error('Error: e_shape must be 3 or 4.');
end
n_np = msh.nbNod;
x_coor = msh.POS(:,1);
y_coor = msh.POS(:,2);

%ID array
ID = zeros(n_np,2);

LM = ID(IEN);





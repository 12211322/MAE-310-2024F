clear
clc
close all
%physical property
v = 0.3;
E = 1e9;
G = E/2/(1+v);

%boundary coundary
fx = @(x,y) 0; fy = @(x,y) 0;
find_xystress; %get the stress function in x_y_coordinate: s_xx, s_yy, s_xy

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
% x_coor = msh.POS(:,1);
% y_coor = msh.POS(:,2);

%true coordinates
pos(:,1) = msh.POS(:,1) + 1;
pos(:,2) = msh.POS(:,2) + 1;

x_coor = pos(:,1);
y_coor = pos(:,2);

%IEN array
IEN_tri = zeros(1,1);
IEN = msh.QUADS(:,1:4);
for ee = 1:size(IEN,1)
    IEN_tri(ee*2-1,1) = IEN(ee,1);
    IEN_tri(ee*2-1,2) = IEN(ee,2);
    IEN_tri(ee*2-1,3) = IEN(ee,3);
    IEN_tri(ee*2,1) = IEN(ee,1);
    IEN_tri(ee*2,2) = IEN(ee,3);
    IEN_tri(ee*2,3) = IEN(ee,4);
end

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

%correct the order of element node 
for i = 1 : n_el/2
    a1 = IEN(i,1); a2 = IEN(i,2); a3 = IEN(i,3); a4 = IEN(i,4); 
    IEN(i,1) = a4; IEN(i,2) = a3; IEN(i,3) = a2; IEN(i,4) = a1; 
    j = i + n_el/2;
    b1 = IEN(j,1); b2 = IEN(j,2); b3 = IEN(j,3); b4 = IEN(j,4); 
    IEN(j,1) = b3; IEN(j,2) = b4; IEN(j,3) = b1; IEN(j,4) = b2; 
    a = IEN(i, 1);
end

%LN array (add the normal vector)
LN = msh.LINES;
for ii = 1: size(LN,1)
    if LN(ii,3) == 10
        LN(ii,4) = 0; LN(ii,5) = 1;
    elseif LN(ii,3) == 8
        LN(ii,4) = 1; LN(ii,5) = 0;
    elseif LN(ii,3) == 12
        coor1 = [pos(LN(ii,1),1),pos(LN(ii,1),2)]; coor2 = [pos(LN(ii,2),1),pos(LN(ii,2),2)];
        mid = (coor1 + coor2)./2;
        LN(ii,4) = -mid(1)/sqrt(mid(1)^2 + mid(2)^2); LN(ii,5) = -mid(2)/sqrt(mid(1)^2 + mid(2)^2);
    elseif LN(ii,3) == 9
        LN(ii,4) = -1; LN(ii,5) = 0;
    elseif LN(ii,3) == 11
        LN(ii,4) = 0; LN(ii,5) = -1;
    end
end


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
def = zeros(n_np,3); %用于储存偏导值

%loop over element to assembly the matrix and vector
for ee = 1:n_el
    x_ele = x_coor(IEN(ee, :));
    y_ele = y_coor(IEN(ee, :));

    k_ele = zeros(2*n_en,2*n_en);
    f_ele = zeros(2*n_en,1);

    for ll = 1:2
        n_x_l = 0.0; n_y_l = 0.0;
        dx_dxi = 0.0; dx_deta = 0.0;
        dy_dxi = 0.0; dy_deta = 0.0;
        xi_n = [-1,1,-1,1]; eta_n = [-1,-1,1,1];
        for aa = 1:n_en %this loop means the shape functions
            n_x_l = n_x_l + x_ele(aa) * Quad(aa, xi_n(ll), eta_n(ll));
            n_y_l = n_y_l + y_ele(aa) * Quad(aa, xi_n(ll), eta_n(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi_n(ll), eta_n(ll));
            dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
            dx_deta = dx_deta + x_ele(aa) * Na_eta;
            dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
            dy_deta = dy_deta + y_ele(aa) * Na_eta;
        end

        % detJ = abs(dx_dxi * dy_deta - dx_deta * dy_dxi); 
        detJ_n = dx_dxi * dy_deta - dx_deta * dy_dxi;%the Jacobian determinant

        for aa = 1:n_en
            n_order = IEN(ee,aa);
            def(n_order,3) = def(n_order,3) + 1;

            Na = Quad(aa, xi_n(ll), eta_n(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi_n(ll), eta_n(ll));
            Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ_n;
            Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ_n;
            
            def(n_order,1) = def(n_order,1) + Na_x;
            def(n_order,2) = def(n_order,2) + Na_y;
        end
    end

    for ll = 1:n_int
        x_l = 0.0; y_l = 0.0; %the coordinate of intergal points
        dx_dxi = 0.0; dx_deta = 0.0;
        dy_dxi = 0.0; dy_deta = 0.0;
        for aa = 1:n_en %this loop means the shape functions
            x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
            y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
            dx_deta = dx_deta + x_ele(aa) * Na_eta;
            dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
            dy_deta = dy_deta + y_ele(aa) * Na_eta;
        end

        % detJ = abs(dx_dxi * dy_deta - dx_deta * dy_dxi); 
        detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;%the Jacobian determinant

        for aa = 1 : n_en
            Na = Quad(aa, xi(ll), eta(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
            Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;

            % f_ele(2*aa - 1) = f_ele(2*aa - 1) + weight(ll) * detJ * fx(x_l, y_l) * Na;
            % f_ele(2*aa) = f_ele(2*aa) + weight(ll) * detJ * fy(x_l, y_l) * Na;
            
            

            for bb = 1 : n_en
                Nb = Quad(bb, xi(ll), eta(ll));
                [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
                Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
                Nb_y = (-Nb_xi * dx_deta + Nb_eta * dx_dxi) / detJ;

                B_a = [Na_x, 0; 0, Na_y; Na_y, Na_x]; B_b = [Nb_x, 0; 0, Nb_y; Nb_y, Nb_x];
                k_ele(2*aa-1:2*aa,2*bb-1:2*bb) = k_ele(2*aa-1:2*aa,2*bb-1:2*bb) + (weight(ll) .* detJ) .* B_a' * D * B_b;


            end % end of bb loop
        end % end of aa loop
    end % end of quadrature loop

    for ii = 1:2
        for aa = 1 : n_en
            pp = 2*(aa-1) + ii;
            PP = LM(ee, pp);
            if PP > 0
                % F(PP) = F(PP) + f_ele(pp);



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

%add the line intergal part on F vector
n_int_bc = 5;
n_en_bc = 2; %number of element nodes
[xi_bc, weight_bc] = Gauss(n_int_bc, -1, 1);
for ee = 1:size(LN,1)
    % if LN(ee,3) == 12
    % else
    f_ele_bc = zeros(2*n_en_bc,1);
    x_ele_bc = x_coor(LN(ee,1:2));
    y_ele_bc = y_coor(LN(ee,1:2));
    s = sqrt((x_ele_bc(1)-x_ele_bc(2))^2 + (y_ele_bc(1)-y_ele_bc(2))^2);
    for ll = 1:n_int_bc
        x_l_bc = 0.0;
        y_l_bc = 0.0;
        for aa = 1:n_en_bc
            x_l_bc = x_l_bc + x_ele_bc(aa) * PolyShape(1, aa, xi_bc(ll), 0);
            y_l_bc = y_l_bc + y_ele_bc(aa) * PolyShape(1, aa, xi_bc(ll), 0);            
        end
        
       
        %y_l_bc = y_ele_bc(1) + (y_ele_bc(2) - y_ele_bc(1))/2*(xi_bc(ll) + 1);

        se_xx = s_xx(x_l_bc,y_l_bc); se_yy = s_yy(x_l_bc,y_l_bc); se_xy = s_xy(x_l_bc,y_l_bc);
        h_v = [se_xx, se_xy; se_xy, se_yy] * [LN(ee,4);LN(ee,5)];
        hx = h_v(1); hy = h_v(2);

        for aa = 1:n_en_bc
            f_ele_bc(2*aa-1) = f_ele_bc(2*aa-1) + weight_bc(ll) * PolyShape(1, aa, xi_bc(ll), 0) * hx * s / 2;
            f_ele_bc(2*aa) = f_ele_bc(2*aa) + weight_bc(ll) * PolyShape(1, aa, xi_bc(ll), 0) * hy * s / 2;
        end
    end

    for dd = 1:2 %degree of freedom
        PP1 = ID(LN(ee,1),dd);
        if PP1 > 0
        F(PP1) = F(PP1) + f_ele_bc(dd);
        end
        PP2 = ID(LN(ee,2),dd);
        if PP2 > 0
        F(PP2) = F(PP2) + f_ele_bc(2+dd);
        end
    end

    
    % end
end


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
    % else
    %     modify disp with the g data. Here it does nothing because g is zero
    % end
end


% strain = zeros(n_el, 3);
% for ee = 1:n_el
% x_ele = [x_coor( IEN(ee, 1:n_en) );x_coor( IEN(ee, 1) )];
% y_ele = [y_coor( IEN(ee, 1:n_en) );y_coor( IEN(ee, 1) )];
% for ll = 1:n_int
% x_l = 0.0; y_l = 0.0;
% dx_dxi = 0.0; dx_deta = 0.0;
% dy_dxi = 0.0; dy_deta = 0.0;
% for aa = 1:n_en
% x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
% y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));
% [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
% dx_dxi = dx_dxi + x_ele(aa) * Na_xi;
% dx_deta = dx_deta + x_ele(aa) * Na_eta;
% dy_dxi = dy_dxi + y_ele(aa) * Na_xi;
% dy_deta = dy_deta + y_ele(aa) * Na_eta;
% end
% detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;
% epsilon = zeros(3, 1);
% for aa = 1:n_en
% Na = Quad(aa, xi(ll), eta(ll));
% [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
% Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
% Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
% B_a = zeros(3, 2);
% B_a(1, 1) = Na_x;
% B_a(2, 2) = Na_y;
% B_a(3, 1) = Na_y;
% B_a(3, 2) = Na_x;
% for i = 1:2
% p = 2 * (aa - 1)+i;
% if i == 1
%     node_disp = x_disp(IEN(ee, aa));
% else
%     node_disp = y_disp(IEN(ee, aa));
% end
% epsilon = epsilon + B_a(:, i) * node_disp;
% end
% end
% % 这里简单取高斯积分点的平均应变作为单元应变
% for qq = 1:3
% strain(ee, qq) = strain(ee, qq) + weight(ll) * epsilon(qq);
% end
% end
% strain(ee, :) = strain(ee, :) / n_int;
% end
% stress = zeros(n_el, 3);
% for ee = 1:n_el
% stress(ee, :) = D * strain(ee, :)';
% end
% % 初始化节点应力数组，每个节点有3个应力分量（xx, yy, xy）
% node_stress = zeros(n_np, 3);
% % 计算每个单元的面积（假设单元为四边形，这里简单用叉积法估算）
% element_area = zeros(n_el, 1);
% for ee = 1:n_el
% x1 = x_coor(IEN(ee, 1));
% y1 = y_coor(IEN(ee, 1));
% x2 = x_coor(IEN(ee, 2));
% y2 = y_coor(IEN(ee, 2));
% x3 = x_coor(IEN(ee, 3));
% y3 = y_coor(IEN(ee, 3));
% x4 = x_coor(IEN(ee, 4));
% y4 = y_coor(IEN(ee, 4));
% area1 = 0.5 * abs((x1 - x3) * (y2 - y4)-(x2 - x4) * (y1 - y3));
% element_area(ee) = area1;
% end
% for ee = 1:n_el
% for aa = 1:n_en
% node_index = IEN(ee, aa);
% % 用单元面积加权
% node_stress(node_index, :) = node_stress(node_index, :) + element_area(ee) * stress(ee, :);
% end
% end
% % 对每个节点的应力进行平均
% for nn = 1:n_np
% total_area = 0;
% for ee = 1:n_el
% if any(IEN(ee, :) == nn)
% total_area = total_area + element_area(ee);
% end
% end
% if total_area > 0
% node_stress(nn, :) = node_stress(nn, :) / total_area;
% end
% end
% 

% exa = zeros(1,1);
% for ii = 1:size(x_disp(:,1))
%     exa(ii) = s_xx(x_coor(ii),y_coor(ii));
% end

hold on
%trisurf(IEN_tri, x_coor, y_coor, node_stress(:, 2));
%trisurf(IEN_tri, x_coor, y_coor, exa);

% axis equal;
% colormap jet
% shading interp
% title('x - direction stress (\sigma_{xx})');
% xlabel('x - coordinate');
% ylabel('y - coordinate');





hold on;

trisurf(IEN_tri, x_coor, y_coor, y_disp);
shading interp;
axis equal;
colormap jet;

























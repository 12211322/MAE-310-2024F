%E = 1e9; v = 0.3; %set the physical property

syms x y E v
xx = x*(1-x)*y*(1-y);
yy = x*(1-x)*y*(1-y);

xx_x = diff(xx,x); xx_y = diff(xx,y);
yy_x = diff(yy,x); yy_y = diff(yy,y);

ep11 = xx_x; ep22 = yy_y; ep12 = (xx_y + yy_x)/2;

%according to plane stress
v_mat = [1, v, 0; v, 1, 0; 0, 0, (1-v)/2];
ep_vect = [ep11; ep22; 2*ep12];
sig_vect = (E/(1-v*v)).*v_mat*ep_vect;
sig11 = sig_vect(1); sig22 = sig_vect(2); sig12 = sig_vect(3);
fx = - diff(sig11,x) - diff(sig12,y);
fy = - diff(sig12,x) - diff(sig22,y);

%according to plane strain
v_mat2 = [v-1, -v, 0; -v, v-1, 0; 0, 0, v-0.5];
sig_vect2 = (E/(v+1)/(2*v-1)).*v_mat2*ep_vect;
sig112 = sig_vect2(1); sig222 = sig_vect2(2); sig122 = sig_vect2(3);
fx2 = - diff(sig112,x) - diff(sig122,y);
fy2 = - diff(sig122,x) - diff(sig222,y);


f_num = subs(fx, v, 0.3);
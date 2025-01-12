syms x y 
T = 1e4;
R = 0.3;
v = 0.3;
E = 1e9;
%D = (E/(1-v*v)).*[1, v, 0; v, 1, 0; 0, 0, 0.5-0.5*v];

r = sqrt(x^2 + y^2);
cosx = x/sqrt(x^2 + y^2);
sinx = y/sqrt(x^2 + y^2);

sig_rr = T/2*(1-(R/r)^2)+T/2*(1-4*((R/r)^2)+3*((R/r)^4))*(cosx^2-sinx^2);
sig_ii = T/2*(1+(R/r)^2)-T/2*(1+3*((R/r)^4))*(cosx^2-sinx^2);
sig_ri = -T/2*(1+2*((R/r)^2)-3*((R/r)^4))*2*sinx*cosx;

sig_xx = sig_rr*(cosx^2) + sig_ii*(sinx^2) - 2*sig_ri*sinx*cosx;
sig_yy = sig_rr*(sinx^2) + sig_ii*(cosx^2) + 2*sig_ri*sinx*cosx;
sig_xy = sig_rr*sinx*cosx - sig_ii*sinx*cosx + sig_ri*(cosx^2 - sinx^2);

sim_xx = simplify(sig_xx);
sim_yy = simplify(sig_yy);
sim_xy = simplify(sig_xy);

% sim_xx = 10000;
% sim_yy = 0;
% sim_xy = 0;

strain = D \ [sim_xx; sim_yy; sim_xy];

u_1 = int(strain(1),'x');
u_2 = int(strain(3)/2 - diff(u_1,'y'),'y');
u_all = u_1 + u_2; %without the integral constant
v_1 = int(strain(2),'y');
v_2 = int(strain(3)/2 - diff(v_1,'x'),'x');
v_all = v_1 + v_2; %without the integral constant

%stress
s_xx = matlabFunction(sim_xx, 'vars', [x, y]);
s_yy = matlabFunction(sim_yy, 'vars', [x, y]);
s_xy = matlabFunction(sim_xy, 'vars', [x, y]);
%strain
strain_xx = matlabFunction(strain(1),'vars',[x, y]);
strain_yy = matlabFunction(strain(2),'vars',[x, y]);
strain_xy = matlabFunction(strain(3)/2,'vars',[x, y]);
%displacement
u_al = matlabFunction(u_all, 'vars',[x, y]);
v_al = matlabFunction(v_all, 'vars',[x, y]);



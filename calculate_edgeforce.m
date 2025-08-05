function [gamma_x,gamma_y,epsilon,gamma1,theta,theta_cos,N]=calculate_edgeforce(eps_xx,eps_yy,eps_xy,Edge_main,v,E,Hx,Hy,d)
 d = d*10^-9;
% eps_xx=eps_xx.*Edge_main;
% eps_xy=eps_xy.*Edge_main;
% eps_yy=eps_yy.*Edge_main;
num = 1;
w=sqrt(Hx.^2+Hy.^2+1);
w_hat=sqrt(Hx.^2+Hy.^2);
w_hat(w_hat==0) = num;
theta_cos = (1-1./w).*Edge_main;%1-cos(theta)
theta = acosd((1./w)).*Edge_main;

N_x = (1./w_hat).*(E/(1-v^2)*Hx.*(eps_xx+v*eps_yy)+E/(1+v)*Hy.*eps_xy);%x_应力
N_y = (1./w_hat).*(E/(1-v^2)*Hy.*(eps_yy+v*eps_xx)+E/(1+v)*Hx.*eps_xy);%y_应力

gamma_x=(1-1./w).*(1./w_hat).*(E/(1-v^2)*Hx.*(eps_xx+v*eps_yy)+E/(1+v)*Hy.*eps_xy);
gamma_y=(1-1./w).*(1./w_hat).*(E/(1-v^2)*Hy.*(eps_yy+v*eps_xx)+E/(1+v)*Hx.*eps_xy);


si = Hy./w_hat;%平面上向量的sin与cos
co = Hx./w_hat;


N = (co.*N_x+si.*N_y)*d.*Edge_main;%正应力乘上了厚度
epsilon = eps_xx.*co.^2+2*si.*co.*eps_xy+si.^2.*eps_yy;%正应变
epsilon = epsilon.*Edge_main;
gamma1 = N.*(1-1./w).*Edge_main;%表观黏附能=正应力*（1-cos(theta)）

gamma_x = gamma_x.*Edge_main;
gamma_y = gamma_y.*Edge_main;
end
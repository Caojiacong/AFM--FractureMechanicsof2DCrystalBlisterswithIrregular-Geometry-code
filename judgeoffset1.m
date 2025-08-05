function [flag]=judgeoffset1(B,N_stress,eps,the0,p,Edge_main,offset,c,theta_c,E,d)
%用于判断矩阵
%change to nm  B:Nm->10^9Nnm  N:N/m->N/nm10^-9
N_stress(N_stress==0) = 1;
flag = Edge_main;
flag(N_stress<0) = 0;
flag(eps<0) = 0;
% N_stress = abs(N_stress);
p= abs(p);
%for the ineq right 
%ineqleft<<offset<<ineqright
theta = the0*pi/180;
ineqleft1 = theta.*(B/E*d).^(1/2).*Edge_main*10^9 - offset/c*Edge_main;
ineqleft = theta.*(B./N_stress).^(1/2).*Edge_main*10^9 - offset/c*Edge_main;%nm%

ineqright = theta.^2.*N_stress./p.*Edge_main*10^9-c*offset.*Edge_main;%nm

flag(ineqleft1>0)=0;
flag(ineqright<0) = 0;
flag(the0<theta_c)=0;
end



 

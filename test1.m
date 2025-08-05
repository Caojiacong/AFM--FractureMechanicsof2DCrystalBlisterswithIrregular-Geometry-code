%关于判断准则的使用结果尚可接受，明确三个步骤，1，处理原始数据，拉平，截取。
% 2，提取出可计算的数据后放入程序进行计算。
% 3，保存计算的数据，利用计算的数据绘制图像，得出结果进行分析
%计算应变时只与面外高度h以及泊松比mu,和远场张力有关
%杨氏模量会影响应力，进而影响切面拉力N
%厚度与弯曲刚度影响计算的压力分布进而影响判断筛选条件，以及远场应力对应变有轻微影响，对于无褶皱样品，我们认为远端张力为0
%%
clear
load("ColorMap_me.mat");
E = 1e12;       %杨氏模量
% E = 0.77e12;
d = 0.335*4;      %厚度
% d = 0.65*4;
poison= 0.168;     %单层石墨烯泊松比
% poison = 0.27;
% B = 2.31e-19;%Nm 单层石墨烯薄膜弯曲刚度
B = E*d^3/(12*(1-poison^2))*10^-27;%Nm 多层石墨烯薄膜弯曲刚度
pathdata = "D:\matlab research\mywork2edtion_improve\EquilibriumofBlisters data_jianguoyun\GronhBN\N=4\sample6round\10291313_540nmround.txt";
% theta1 = 28.1*pi/180;
% theta2 = -14*pi/180;
% theta3 = 75.1*pi/180;
% theta1 = 11.5*pi/180;
% theta2 = -29.2*pi/180;
% theta3 = -94.1*pi/180;
% theta1 = -18.4*pi/180;
% theta2 = 18.8*pi/180;
% theta3 = -49.6*pi/180;
len = 540;
%%
% [Ksi] = get_forcefiled(line1,line2,line3,theta1,theta2,theta3,B,d,E);
[Ksi] = [0 0 0 ];
% %%
FAMA=load(pathdata);%加载数据
FAMA=FAMA/(1e-9); %转换为nm
N=120;%Chebyshev点数目
Numofgrid = 800; %计算划分网格数
sigma = 13;     %平滑参数 gr/gr13 gr/hBN19
c_x=Ksi(1);%Tx/E
c_y=Ksi(3);%Ty/E
c_xy=Ksi(2);%Txy/E
% c_x =-1.577E-5; c_y = 1.646E-7; c_xy= 4.1345E-6;
[imu_f,imL_u]=dataprocessing_new(FAMA,len,len,sigma,Numofgrid);


[eps_xx,eps_yy,eps_xy,chi,gx,gy,~]=straintensor_solver(imu_f,len,len,N,poison,c_x,c_y,c_xy);
plotresult(imL_u,imu_f,chi,eps_xx,eps_yy,eps_xy,len,len,ColorMap_me);
 %%
 %计算面外应力
[p,Hx,Hy,doublelaplace_h,P,p1,p2]=solve_P(imu_f,len,d,poison,chi,E);
px=-p.*Hx;
py=-p.*Hy;
%%
h_max = max(imu_f(:));
offset=0.15*h_max;%临界高度
%%求给定形貌图的气泡边界



[flag1,flag_main,Edge_main]=regionalism(imu_f,offset);
[row,col] = find(Edge_main);
 r_a = max(sqrt((row-Numofgrid/2).^2+(col-Numofgrid/2).^2))*len/(Numofgrid-1);
 P_app = E*d*h_max^3/r_a^4;

[gamma_x,gamma_y,eps,gamma,the0,theta_cos,N_stress]=calculate_edgeforce(eps_xx,eps_yy,eps_xy,Edge_main,poison,E,Hx,Hy,d);
clear row col line1 line2 lin3 theta1 theta2 theta3 H 
pmax = 1/4*max(p(:));
judgeflag = judgeoffset1(B,N_stress,eps,the0,pmax,Edge_main,offset,5,8,E,d);
scatter(the0(judgeflag==1),gamma(judgeflag==1))
%scatter(the0(the0~=0),gamma(the0~=0))
% save("D:\matlab research\mywork2edtion_improve\code\plot result\dataset\Gr on hBN\4L\4GronhBNround6.mat")
function [imu_f,imL_u]=dataprocessing_new(FAMA,xlen,ylen,sigma,Numofgrid)
%未解决矩阵不是方阵的困难
[yd,xd]=size(FAMA);%yd为传入矩阵的行数，xd为传入矩阵的列数
dx0=xlen/(xd-1);dy0=ylen/(yd-1);%分别为x方向与y方向网格间距
x=dx0*(0:xd-1)';
y=dy0*(0:yd-1)';
x=x-median(x);
y=y-median(y);%将求解的物理区域的中心作为坐标原点。
[X0,Y0]=meshgrid(x,y);%给定物理区域的节点
xu=linspace(min(x),max(x),Numofgrid);
yu=linspace(min(y),max(y),Numofgrid);
[Xu,Yu]=meshgrid(xu,yu);
%样条插值
imL_u=interp2(X0,Y0,FAMA,Xu,Yu,'spline');
H = max(imL_u(:));
% 
%%%光滑处理 
imu_f=imgaussfilt(imL_u,sigma);
% [imu_f] = smoothm(imL_u,sigma,xlen);
imu_f = imu_f/max(imu_f(:))*H;
%,"FilterSize",5
end
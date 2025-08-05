function [eps_xx,eps_yy,eps_xy,chi,gx,gy,H]=straintensor_solver(imu_f,xlen,ylen,N,poison,c_x,c_y,c_xy)
%where c_x=sigma_x/E  c_y=sigma_y/E  c_xy=gama_xy/E
%% Gauss-Lobbatto Points and Chebyshev derivative matrix 
[D,gl] = chebyshev(N);%求得chebyshev谱求导矩阵
D2 = D^2;
% Define 2D Laplacian Operator Using kron.m take compute the tensor product
% of the second derivative matrix with the identiy matrix. This creates the
% Laplcian which is an (N+1)^2 x (N+1)^2 matrix for the whole space.
Id = eye(N+1);
L = (sparse(kron(Id,D2)+kron(D2,Id)));%拉普拉斯算符

%% Source function (gaussian curvature)
[Xgl,Ygl] = meshgrid(gl,gl);

xx=Xgl(:); yy = Ygl(:); % colon operator here stacks all columns into a 1D vector.

H = 2*imu_f/xlen;%无量纲化    %%%&&&   注意这在无量纲化时若x，y长度不同要做更改

[gx,gy] = gradient(H,2/(size(H,1)-1));%%
[gxx,gxy] = gradient(gx,2/(size(H,1)-1));%%
[~,gyy] = gradient(gy,2/(size(H,1)-1));%%

xi = linspace(-1,1,size(H,1));
yi = linspace(-1,1,size(H,1));
[Xi,Yi] = meshgrid(xi,yi);

% Use 2D interpolation to calculate the curvature at the collocation points
gxxi = interp2(Xi,Yi,gxx,Xgl,Ygl); 
gyyi = interp2(Xi,Yi,gyy,Xgl,Ygl);
gxyi = interp2(Xi,Yi,gxy,Xgl,Ygl);

f = gxxi.*gyyi-(gxyi).^2;
f = -f(:);
%修改谱求导矩阵

bound = find(abs(xx)==1 | abs(yy)==1);
L(bound,:) = 0;L(bound,bound) = eye(length(bound));
f(bound) = c_x+c_y; %%%     

tic

v = L\f;
v(bound) = 1/2*c_x*yy(bound).^2+1/2*c_y*xx(bound).^2 ...
        -c_xy*xx(bound).*yy(bound);%%%
Phi = v;
chi_t = L\v;
toc
chi_t = gather(chi_t);
Ni = size(H,1);%对求出来的应力函数插值，插值点数目为Ni即为 变量xu的数目
xf = linspace(-1,1,Ni);
yf = linspace(-1,1,Ni);
[Xf,Yf] = meshgrid(xf,yf);

chi = reshape(chi_t,[N+1,N+1]);
chi = interp2(Xgl,Ygl,chi,Xf,Yf,'spline');%4x/(l^2*E)
dx = 2/(Ni-1);
[chi_x,chi_y] = gradient(chi,dx);
[chi_xx,chi_xy] = gradient(chi_x,dx);
[~,chi_yy] = gradient(chi_y,dx);

eps_xx = chi_yy-poison*chi_xx;
eps_yy = chi_xx-poison*chi_yy;
eps_xy = -(1+poison)*chi_xy;%%



end



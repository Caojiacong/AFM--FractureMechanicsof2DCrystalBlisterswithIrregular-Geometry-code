function plotresult(imL_u,imu_f,chi,eps_xx,eps_yy,eps_xy,xlen,ylen,ColorMap_me)

[Nx,Ny]=size(imu_f);
xf=linspace(-1,1,Nx);
yf=linspace(-1,1,Ny);
[Xf,Yf]=meshgrid(xf,yf);
figure(1)
imagesc(xf*xlen,yf*ylen,chi), axis image ,colorbar('Southoutside')
title('归一化 Airy stress function')

figure(2), 
set(gcf,'position',[250 300 1000 400])
p1=subplot(131);
imagesc(xf*xlen/2,xf*ylen/2,imu_f), axis image, colorbar('eastoutside',FontSize=13)
set(p1,'FontSize',13)
title('Nanobubble Topography','FontSize',16)
p2=subplot(132); 
imagesc(xf*xlen/2,xf*ylen/2,(eps_xx + eps_yy)), axis image, colorbar('eastoutside',FontSize=13), % multiply by 100 to convert to %
set(p2,'FontSize',13)
title('\epsilon_{\fontsize{15}\itxx} + \epsilon_{\fontsize{15}\ityy}','FontSize',20) 
p3=subplot(133);imagesc(xf*xlen/2,xf*xlen/2,eps_xy),axis image, colorbar('eastoutside',FontSize=16)
set(p3,'FontSize',13)
title('\fontsize{20}\epsilon_{\fontsize{15}\itxy}')
colormap(p1,'hot')
colormap(p2,ColorMap_me)
colormap(p3,ColorMap_me)
figure(4)
surf(Xf*xlen/2,Yf*ylen/2,imu_f,'EdgeColor','interp');
% surf(Xf*xlen/2,Yf*ylen/2,imu_f);
title('形貌图')
xlabel('X/nm') , ylabel('Y/nm')

figure(5)
surf(Xf*xlen/2,Yf*ylen/2,(eps_xx+eps_yy),'EdgeColor','interp'); title('trace 高度图')

figure(6)

p4 = subplot(121);
imagesc(xf*xlen/2,yf*ylen/2,imL_u)
set(p4,'FontSize',18)
title('Before Smooth',FontSize=20)
axis image, colormap hot; colorbar()
% axis off
 colorbar;
c = colorbar;
set(c,'Ticks',0:floor(1/5*max(imu_f(:))):max(imu_f(:)),'FontSize',18) 
p5 = subplot(122);
imagesc(xf*xlen/2,yf*ylen/2,imu_f)
set(p5,'FontSize',18)
title('After Smooth','FontSize',20)
axis image, colormap hot; 
colorbar;
c = colorbar;
set(c,'Ticks',0:floor(1/5*max(imu_f(:))):max(imu_f(:)),'FontSize',18) 

figure(59)
surf(Xf*xlen/2,Yf*ylen/2,imL_u,'EdgeColor','interp');
% surf(Xf*xlen/2,Yf*ylen/2,imu_f);
title('未处理形貌图')
xlabel('X/nm') , ylabel('Y/nm')
end
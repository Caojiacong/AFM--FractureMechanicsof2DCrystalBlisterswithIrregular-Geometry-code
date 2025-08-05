function [flag1,flag_main,Edge_main]=regionalism(imu_f,offset)
%%
%%求给定形貌图的主要大气泡的边界
[Nx,Ny]=size(imu_f);
% dx=len/Ny;
flag=zeros(Nx,Ny);

for i=1:Nx
    for j=1:Ny
        if imu_f(i,j)>offset
            flag(i,j)=1;
        end
    end
end

[flag1,NUM] = bwlabel(flag,8);%寻找连通区域
numofblock = zeros(1,NUM);
%寻找最大连通区域，maxblock对应于flag1中最大连通区域的编号，
for k = 1:NUM
    count = 0;
    for i = 1:Nx
        for j = 1:Ny
            if flag1(i,j) == k
                count = count+1;
            end
        end
    end
    numofblock(1,k) = count;
end
maxblock = 1;
for i = 1:NUM
    if numofblock(i) > numofblock(maxblock)
        maxblock = i;
    end
end

flag_main =zeros(Nx,Ny);
% flag2=zeros(Nx,Ny);
% uflag=zeros(Nx,Ny);
% dflag=zeros(Nx,Ny);

for i=1:Nx
    for j=1:Ny
        if flag1(i,j) == maxblock
            
            flag_main(i,j) = 1;
           
        end
        
    end
end

Edge_main=edge(flag_main);
clear i j k
end
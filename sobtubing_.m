function sobtubing_()
tic;close all
ee=1e-8;
%划分时空网格
dx=0.01;dy=0.01; dt=0.0002;
x=0:dx:4;
y=0:dy:1;
Nx=length(x);
Ny=length(y);%定义空间网格数目
%一维sod激波管时空网格划分
%A=2.1916;
%while A*dt/dx>=0.5
%    dt=dt/10;
%end
Nt=0.2/dt;
%Nright=round(1/dx+1);Nleft=round(9/dx+1);
%disp(Nright);
%disp(Nleft);
%N=Nright+Nleft-1;
%初始条件
%P=[10.333.*ones(1,Nleft-1) ones(1,Nright)];
%Den=[3.587.*ones(1,Nleft-1) ones(1,Nright)+0.2.*sin(5.*x1)];
%u=[zeros(1,Nleft-1) zeros(1,Nright)];
%v=[zeros(1,Nleft-1) zeros(1,Nright)];
P=ones(Ny,Nx);
Den=ones(Ny,Nx);
u=ones(Ny,Nx);
v=ones(Ny,Nx);
% %二维sod问题的初始值
% for i=1:201
%     for j=1:201
%      if sqrt(x(1,i)^2+y(1,j)^2)<=0.4
%         P(i,j)=1;
%         Den(i,j)=1;
%         u(i,j)=0;
%         v(i,j)=0;
%      else
%         P(i,j)=0.1;
%         Den(i,j)=0.125;
%         u(i,j)=0;
%         v(i,j)=0;
%      end
%     end
% end
%双马赫反射问题初始值
for m=1:Ny
    for n=1:Nx
        if sqrt(3)*(x(1,n)-1/6)-y(1,m)<=0
             Den(m,n)=8;u(m,n)=4.125*sqrt(3);v(m,n)=-4.125;P(m,n)=116.5;%赋予初始             
        else 
             Den(m,n)=1.4;u(m,n)=0;v(m,n)=0;P(m,n)=1;
        end
    end
end
gama=1.4;
Den_u=Den.*u;
Den_v=Den.*v;
E=P./(gama-1)+(0.5.*(u.^2+v.^2)).*Den;
%[~,s]=contourf(x,y,Den);
%s.LineColor = 'none';
%shading interp
%axis equal
%title("密度分布图")
%xlabel("X")
%ylabel("Y")
%初始化边界条件
%左边界
%for k=1:Ny
%    Den(k,1)=8;u(k,1)=4.125*sqrt(3);v(k,1)=-4.125;P(k,1)=116.5;
%end
%右边界
%for k=1:Ny
%    Den(k,Nx)=1.4;u(k,Nx)=0;v(k,Nx)=0;P(k,Nx)=1;
%end  
%下边界
%for j=2:Nx-1
% if x(1,j)<=1/6
%     Den(1,j)=8;u(1,j)=4.125*sqrt(3);v(1,j)=-4.125;P(1,j)=116.5;
% else
%     Den(1,j)=1.4;u(1,j)=0;v(1,j)=0;P(1,j)=1;
% end
%end
%开始计算流场
for t=1:Nt
    %进入时间步，进行变量初始化
    Den_u=Den.*u;
    Den_v=Den.*v;
    E=P./(gama-1)+(0.5.*(u.^2+v.^2)).*Den;
    %进行分裂变量计算
    epso=1e-8*ones(Ny,Nx);
    gama=1.4;
    C=sqrt(gama*P./Den);
    lamtax1=u;lamtax2=u-C;lamtax3=u+C;
    lamtay1=v;lamtay2=v-C;lamtay3=v+C;
    lamta_px1=0.5.*(lamtax1+sqrt(lamtax1.^2+epso.^2));
    lamta_nx1=0.5.*(lamtax1-sqrt(lamtax1.^2+epso.^2));
    lamta_py1=0.5.*(lamtay1+sqrt(lamtay1.^2+epso.^2));
    lamta_ny1=0.5.*(lamtay1-sqrt(lamtay1.^2+epso.^2));
    lamta_px2=0.5.*(lamtax2+sqrt(lamtax2.^2+epso.^2));
    lamta_nx2=0.5.*(lamtax2-sqrt(lamtax2.^2+epso.^2));
    lamta_py2=0.5.*(lamtay2+sqrt(lamtay2.^2+epso.^2));
    lamta_ny2=0.5.*(lamtay2-sqrt(lamtay2.^2+epso.^2));
    lamta_px3=0.5.*(lamtax3+sqrt(lamtax3.^2+epso.^2));
    lamta_nx3=0.5.*(lamtax3-sqrt(lamtax3.^2+epso.^2));
    lamta_py3=0.5.*(lamtay3+sqrt(lamtay3.^2+epso.^2));
    lamta_ny3=0.5.*(lamtay3-sqrt(lamtay3.^2+epso.^2));
    %计算正通量
    gama=1.4;
    C=sqrt(gama.*P./Den);
    f1_Posx=0.5/gama.*Den.*(2*(gama-1).*lamta_px1+lamta_px2+lamta_px3);
    f2_Posx=0.5/gama.*Den.*(2*(gama-1).*lamta_px1.*u+lamta_px2.*(u-C)+lamta_px3.*(u+C));
    f3_Posx=0.5/gama.*Den.*(2*(gama-1).*lamta_px1.*v+lamta_px2.*v+lamta_px3.*v);
    f4_Posx=0.5/gama.*Den.*((gama-1).*lamta_px1.*(u.^2+v.^2)+0.5.*lamta_px2.*((u-C).^2+v.^2)+0.5.*lamta_px3.*((u+C).^2+v.^2)+(0.5.*(3-gama)/(gama-1).*(lamta_px2+lamta_px3).*C.^2));
    f1_Posy=0.5/gama.*Den.*(2*(gama-1).*lamta_py1+lamta_py2+lamta_py3);
    f2_Posy=0.5/gama.*Den.*(2*(gama-1).*lamta_py1.*u+lamta_py2.*u+lamta_py3.*u);
    f3_Posy=0.5/gama.*Den.*(2*(gama-1).*lamta_py1.*v+lamta_py2.*(v-C)+lamta_py3.*(v+C));
    f4_Posy=0.5/gama.*Den.*((gama-1).*lamta_py1.*(u.^2+v.^2)+0.5*lamta_py2.*((v-C).^2+u.^2)+0.5*lamta_py3.*((v+C).^2+u.^2)+(0.5*(3-gama)/(gama-1).*(lamta_py2+lamta_py3).*C.^2));
    %计算负通量
    gama=1.4;
    C=sqrt(gama*P./Den);
    f1_Negx=0.5/gama.*Den.*(2*(gama-1).*lamta_nx1+lamta_nx2+lamta_nx3);
    f2_Negx=0.5/gama.*Den.*(2*(gama-1).*lamta_nx1.*u+lamta_nx2.*(u-C)+lamta_nx3.*(u+C));
    f3_Negx=0.5/gama.*Den.*(2*(gama-1).*lamta_nx1.*v+lamta_nx2.*v+lamta_nx3.*v);
    f4_Negx=0.5/gama.*Den.*((gama-1).*lamta_nx1.*(u.^2+v.^2)+0.5.*lamta_nx2.*((u-C).^2+v.^2)+0.5.*lamta_nx3.*((u+C).^2+v.^2)+(0.5.*(3-gama)/(gama-1).*(lamta_nx2+lamta_nx3).*C.^2));
    f1_Negy=0.5/gama.*Den.*(2*(gama-1).*lamta_ny1+lamta_ny2+lamta_ny3);
    f2_Negy=0.5/gama.*Den.*(2*(gama-1).*lamta_ny1.*u+lamta_ny2.*u+lamta_ny3.*u);
    f3_Negy=0.5/gama.*Den.*(2*(gama-1).*lamta_ny1.*v+lamta_ny2.*(v-C)+lamta_ny3.*(v+C));
    f4_Negy=0.5/gama.*Den.*((gama-1).*lamta_ny1.*(u.^2+v.^2)+0.5.*lamta_ny2.*((v-C).^2+u.^2)+0.5*lamta_ny3.*((v+C).^2+u.^2)+(0.5*(3-gama)/(gama-1).*(lamta_ny2+lamta_ny3).*C.^2));
    %内部流场计算      
    for i=2:Ny-1
     for j=2:Nx-1
      %计算密度
      temp1x=(f1_Posx(i,j)-f1_Posx(i,j-1))/dx;
      temp2x=(f1_Negx(i,j+1)-f1_Negx(i,j))/dx;
      temp1y=(f1_Posy(i,j)-f1_Posy(i-1,j))/dy;
      temp2y=(f1_Negy(i+1,j)-f1_Negy(i,j))/dy;
      Den(i,j)=Den(i,j)-dt*(temp1x+temp2x+temp1y+temp2y);
      %密度速度x乘积计算
      temp1x=(f2_Posx(i,j)-f2_Posx(i,j-1))/dx;
      temp2x=(f2_Negx(i,j+1)-f2_Negx(i,j))/dx;
      temp1y=(f2_Posy(i,j)-f2_Posy(i-1,j))/dy;
      temp2y=(f2_Negy(i+1,j)-f2_Negy(i,j))/dy;
      Den_u(i,j)=Den_u(i,j)-dt*(temp1x+temp2x+temp1y+temp2y);
      %计算速度x
      u(i,j)=Den_u(i,j)/Den(i,j);
      %密度速度y乘积计算
      temp1x=(f3_Posx(i,j)-f3_Posx(i,j-1))/dx;
      temp2x=(f3_Negx(i,j+1)-f3_Negx(i,j))/dx;
      temp1y=(f3_Posy(i,j)-f3_Posy(i-1,j))/dy;
      temp2y=(f3_Negy(i+1,j)-f3_Negy(i,j))/dy;
      Den_v(i,j)=Den_v(i,j)-dt*(temp1x+temp2x+temp1y+temp2y);
      %计算速度y
      v(i,j)=Den_v(i,j)/Den(i,j);
      %能量计算
      temp1x=(f4_Posx(i,j)-f4_Posx(i,j-1))/dx;
      temp2x=(f4_Negx(i,j+1)-f4_Negx(i,j))/dx;
      temp1y=(f4_Posy(i,j)-f4_Posy(i-1,j))/dy;
      temp2y=(f4_Negy(i+1,j)-f4_Negy(i,j))/dy;
      E(i,j)=E(i,j)-dt*(temp2x+temp1x+temp2y+temp1y);
      %压强计算
      P(i,j)=(gama-1)*(E(i,j)-0.5*Den(i,j)*(u(i,j)^2+v(i,j)^2));
     end
    end
    %边界流场计算
    %上边界流场计算
    for j=2:Nx-1
     %计算密度
     temp1x=(f1_Posx(Ny,j)-f1_Posx(Ny,j-1))/dx;
     temp2x=(f1_Negx(Ny,j+1)-f1_Negx(Ny,j))/dx;
     temp1y=(f1_Posy(Ny,j)-f1_Posy(Ny-1,j))/dy;
     Den(Ny,j)=Den(Ny,j)-dt*(temp1x+temp2x+temp1y);
   %密度速度x乘积计算
     temp1x=(f2_Posx(Ny,j)-f2_Posx(Ny,j-1))/dx;
     temp2x=(f2_Negx(Ny,j+1)-f2_Negx(Ny,j))/dx;
     temp1y=(f2_Posy(Ny,j)-f2_Posy(Ny-1,j))/dy;
     Den_u(Ny,j)=Den_u(Ny,j)-dt*(temp1x+temp2x+temp1y);
     %计算速度x
     u(Ny,j)=Den_u(Ny,j)/Den(Ny,j);
     %密度速度y乘积计算
     temp1x=(f3_Posx(Ny,j)-f3_Posx(Ny,j-1))/dx;
     temp2x=(f3_Negx(Ny,j+1)-f3_Negx(Ny,j))/dx;
     temp1y=(f3_Posy(Ny,j)-f3_Posy(Ny-1,j))/dy;
     Den_v(Ny,j)=Den_v(Ny,j)-dt*(temp1x+temp2x+temp1y);
     %计算速度y
     v(Ny,j)=Den_v(Ny,j)/Den(Ny,j);
     %能量计算
     temp1x=(f4_Posx(Ny,j)-f4_Posx(Ny,j-1))/dx;
     temp2x=(f4_Negx(Ny,j+1)-f4_Negx(Ny,j))/dx;
     temp1y=(f4_Posy(Ny,j)-f4_Posy(Ny-1,j))/dy;      
     E(Ny,j)=E(Ny,j)-dt*(temp2x+temp1x+temp1y);
     %压强计算
    P(Ny,j)=(gama-1)*(E(Ny,j)-0.5*Den(Ny,j)*(u(Ny,j)^2+v(Ny,j)^2));
    end
    %右边界流场计算
%    for j=2:Ny-1
%      %计算密度
%      temp1x=(f1_Posx(j,Nx)-f1_Posx(j,Nx-1))/dx;
%      temp1y=(f1_Posy(j,Nx)-f1_Posy(j-1,Nx))/dy;
%      temp2y=(f1_Negy(j+1,Nx)-f1_Negy(j,Nx))/dy;
%      Den(j,Nx)=Den(j,Nx)-dt*(temp1x+temp1y+temp2y);
%      %密度速度x乘积计算
%      temp1x=(f2_Posx(j,Nx)-f2_Posx(j,Nx-1))/dx;
%      temp1y=(f2_Posy(j,Nx)-f2_Posy(j-1,Nx))/dy;
%      temp2y=(f2_Negy(j+1,Nx)-f2_Negy(j,Nx))/dy;
%      Den_u(j,Nx)=Den_u(j,Nx)-dt*(temp1x+temp1y+temp2y);
%      %计算速度x
 %     u(j,Nx)=Den_u(j,Nx)/Den(j,Nx);
 %     %密度速度y乘积计算
%      temp1x=(f3_Posx(j,Nx)-f3_Posx(j,Nx-1))/dx;
%      temp1y=(f3_Posy(j,Nx)-f3_Posy(j-1,Nx))/dy;
%      temp2y=(f3_Negy(j+1,Nx)-f3_Negy(j,Nx))/dy;
%      Den_v(j,Nx)=Den_v(j,Nx)-dt*(temp1x+temp1y+temp2y);
%      %计算速度y
%      v(j,Nx)=Den_v(j,Nx)/Den(j,Nx);
%      %能量计算
%      temp1x=(f4_Posx(j,Nx)-f4_Posx(j,Nx-1))/dx;
%      temp1y=(f4_Posy(j,Nx)-f4_Posy(j-1,Nx))/dy;
%      temp2y=(f4_Negy(j+1,Nx)-f4_Negy(j,Nx))/dy;
%      E(j,Nx)=E(j,Nx)-dt*(temp1x+temp2y+temp1y);
%      %压强计算
%      P(j,Nx)=(gama-1)*(E(j,Nx)-0.5*Den(j,Nx)*(u(j,Nx)^2+v(j,Nx)^2));
%    end
%    %左边界
%    for k=1:Ny
%        Den(k,1)=8;u(k,1)=4.125*sqrt(3);v(k,1)=-4.125;P(k,1)=116.5;
%    end
    %下边界流场计算
    for j=2:Nx-1
      if x(1,j)<1/6
      Den(j,1)=8;u(j,1)=4.125*sqrt(3);v(j,1)=-4.125;P(j,1)=116.5;
      else
      %计算密度
       temp1x=(f1_Posx(1,j)-f1_Posx(1,j-1))/dx;
       temp2x=(f1_Negx(1,j+1)-f1_Negx(1,j))/dx;
       temp2y=(f1_Negy(2,j)-f1_Negy(1,j))/dy;
       Den(1,j)=Den(1,j)-dt*(temp1x+temp2x+temp2y);
       %密度速度x乘积计算
       temp1x=(f2_Posx(1,j)-f2_Posx(1,j-1))/dx;
       temp2x=(f2_Negx(1,j+1)-f2_Negx(1,j))/dx;
       temp2y=(f2_Negy(2,j)-f2_Negy(1,j))/dy;
       Den_u(1,j)=Den_u(1,j)-dt*(temp1x+temp2x+temp2y);
       %计算速度x
       u(1,j)=Den_u(1,j)/Den(1,j);
       %密度速度y乘积计算
       temp1x=(f3_Posx(1,j)-f3_Posx(1,j-1))/dx;
       temp2x=(f3_Negx(1,j+1)-f3_Negx(1,j))/dx;
       temp2y=(f3_Negy(2,j)-f3_Negy(1,j))/dy;
       Den_v(1,j)=Den_v(1,j)-dt*(temp1x+temp2x+temp2y);
       %计算速度y
       v(1,j)=Den_v(1,j)/Den(1,j);
       %能量计算
       temp1x=(f4_Posx(1,j)-f4_Posx(1,j-1))/dx;
       temp2x=(f4_Negx(1,j+1)-f4_Negx(1,j))/dx;
       temp2y=(f4_Negy(2,j)-f4_Negy(1,j))/dy;
       E(1,j)=E(1,j)-dt*(temp2x+temp1x+temp2y);
       %压强计算
       P(1,j)=(gama-1)*(E(1,j)-0.5*Den(1,j)*(u(1,j)^2+v(1,j)^2));
       end
     end
     %反射
     for i=1:Nx
         if x(1,i)>1/6&&v(1,i)<=0
             v(1,i)=-v(1,i);
        end
     end
end 
%绘图
[~,s]=contourf(x,y,P);
s.LineColor = 'none';
shading interp
axis equal
title("压力分布图")
xlabel("X")
ylabel("Y")
Calculate_time=toc;
disp(Calculate_time);
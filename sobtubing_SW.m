function sobtubing_SW()
tic;close all
ee=1e-8;
%划分时空网格
dx=0.01; dt=0.1;
A=2.1916;
while A*dt/dx>=0.5
    dt=dt/10;
end
Nt=round(2/dt);
Nright=round(5/dx+1);Nleft=round(5/dx+1);
N=Nright+Nleft-1;
%初始条件
P=[ones(1,Nleft-1) 0.1*ones(1,Nright)];
Den=[ones(1,Nleft-1) 0.125*ones(1,Nright)];
u=zeros(1,N);
gama=1.4;
Den_u=Den.*u;
E=P./(gama-1)+(0.5*u.^2).*Den;
%计算特征分裂
for j=1:Nt
    epso=1e-8*ones(1,N);
    gama=1.4;
    C=sqrt(gama*P./Den);
    lamta=ones(3,N);lamta_p=ones(3,N);lamta_n=ones(3,N);
    lamta(1,:)=u;lamta(2,:)=u-C;lamta(3,:)=u+C;
    for i=1:3
        lamta_p(i,:)=0.5*(lamta(i,:)+sqrt(lamta(i,:).^2+epso.^2));
        lamta_n(i,:)=0.5*(lamta(i,:)-sqrt(lamta(i,:).^2+epso.^2));
    end
    %计算正通量
    gama=1.4;
    C=sqrt(gama*P./Den);
    Ftran_p=ones(3,N);
    f1_Pos=0.5/gama*Den.*(2*(gama-1)*lamta_p(1,:)+lamta_p(2,:)+lamta_p(3,:));
    f2_Pos=0.5/gama*Den.*(2*(gama-1)*lamta_p(1,:).*u+lamta_p(2,:).*(u-C)+lamta_p(3,:).*(u+C));
    f3_Pos=0.5/gama*Den.*((gama-1)*lamta_p(1,:).*u.^2+0.5*lamta_p(2,:).*(u-C).^2+0.5*lamta_p(3,:).*(u+C).^2+(0.5*(3-gama)/(gama-1)*(lamta_p(2,:)+lamta_p(3,:)).*C.^2));
    %计算负通量
    gama=1.4;
    C=sqrt(gama*P./Den);
    f1_Neg=0.5/gama*Den.*(2*(gama-1)*lamta_n(1,:)+lamta_n(2,:)+lamta_n(3,:));
    f2_Neg=0.5/gama*Den.*(2*(gama-1)*lamta_n(1,:).*u+lamta_n(2,:).*(u-C)+lamta_n(3,:).*(u+C));
    f3_Neg=0.5/gama*Den.*((gama-1)*lamta_n(1,:).*u.^2+0.5*lamta_n(2,:).*(u-C).^2+0.5*lamta_n(3,:).*(u+C).^2+(0.5*(3-gama)/(gama-1)*(lamta_n(2,:)+lamta_n(3,:)).*C.^2));
    %计算流动参数
    for i=2:N-1
     %计算密度
     temp1(1,i)=(f1_Pos(1,i)-f1_Pos(1,i-1))/dx;
     temp2(1,i)=(f1_Neg(1,i+1)-f1_Neg(1,i))/dx;
     Den(1,i)=Den(1,i)-dt*(temp1(1,i)+temp2(1,i));
     %密度速度乘积计算
     temp1(1,i)=(f2_Pos(1,i)-f2_Pos(1,i-1))/dx;
     temp2(1,i)=(f2_Neg(1,i+1)-f2_Neg(1,i))/dx;
     Den_u(1,i)=Den_u(1,i)-dt*(temp2(1,i)+temp1(1,i));
     %计算速度
     u(1,i)=Den_u(1,i)/Den(1,i);
     %能量计算
     temp1(1,i)=(f3_Pos(1,i)-f3_Pos(1,i-1))/dx;
     temp2(1,i)=(f3_Neg(1,i+1)-f3_Neg(1,i))/dx;
     E(1,i)=E(1,i)-dt*(temp2(1,i)+temp1(1,i));
     %压强计算
     P(1,i)=(gama-1)*(E(1,i)-0.5*Den(1,i)*u(1,i)^2);
    end
end
%绘图
x=-5:0.01:5;
axis([0 1 0 1]);
plot(x,u,'Linewidth',1.2,'Color','r');hold on;
plot(x,P,'Linewidth',1.2,'Color','g');hold on;
plot(x,Den,'Linewidth',1.2,'Color','b');hold off;
legend('t=2速度分布','t=2压力分布','t=2密度分布')
Calculate_time=toc
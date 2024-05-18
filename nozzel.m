function SW()
clc;clear;
tic;close all
ee=1e-8;
%划分时空网格（C=0.5）
dx=0.1; dt=0.0207;
x=0:dx:3;
t=0:dt:28.952;
Nx=length(x);
Nt=length(t);
%初始条件（根据文中7-74）
P=ones(1,Nx)-0.314*x.*ones(1,Nx);
T=ones(1,Nx)-0.2314*x.*ones(1,Nx);
A=ones(1,Nx)+2.2.*((x-1.5.*ones(1,Nx)).^2).*ones(1,Nx);%定义喷管7-73
u=(0.1.*ones(1,Nx)+1.09.*x).*sqrt(T);
%定义偏导数以及预估步修正步矩阵
P1=zeros(1,Nx);u1=zeros(1,Nx);T1=zeros(1,Nx);
P2=zeros(1,Nx);u2=zeros(1,Nx);T2=zeros(1,Nx);
P3=zeros(1,Nx);u3=zeros(1,Nx);T3=zeros(1,Nx);
P4=zeros(1,Nx);u4=zeros(1,Nx);T4=zeros(1,Nx);
Ma=u./sqrt(T);%mach number定义
%计算
for j=1:Nt%时间
    gama=1.4;
    P2(1,1)=1;T2(1,1)=1;u2(1,1)=u(1,1);%定义修正步（1，1）位置数值
    for i=2:Nx-1
        %计算向后差分偏导数（7-51，7-52，7-53）%ps:7-52式子有错，应当加负号
        P1(1,i)=-u(1,i)*(P(1,i+1)-P(1,i))/dx-P(1,i)*(u(1,i+1)-u(1,i))/dx-P(1,i)*u(1,i)*(log(A(1,i+1))-log(A(1,i)))/dx;
        u1(1,i)=-u(1,i)*(u(1,i+1)-u(1,i))/dx-(1/gama)*(((T(1,i+1)-T(1,i))/dx)+(T(1,i)*(P(1,i+1)-P(1,i))/(P(1,i)*dx)));
        T1(1,i)=-u(1,i)*(T(1,i+1)-T(1,i))/dx-(gama-1)*T(1,i)*((u(1,i+1)-u(1,i))/dx+u(1,i)*(log(A(1,i+1))-log(A(1,i)))/dx);
        %计算预估值（7-54，7-55，7-56）
        P2(1,i)=P(1,i)+dt*P1(1,i);
        u2(1,i)=u(1,i)+dt*u1(1,i);
        T2(1,i)=T(1,i)+dt*T1(1,i);
        %计算向前差分（7-57，7-58，7-59）
        P3(1,i)=-u2(1,i)*(P2(1,i)-P2(1,i-1))/dx-P2(1,i)*(u2(1,i)-u2(1,i-1))/dx-P2(1,i).*u2(1,i)*(log(A(1,i))-log(A(1,i-1)))/dx;
        u3(1,i)=-u2(1,i)*(u2(1,i)-u2(1,i-1))/dx-(1/gama)*((T2(1,i)-T2(1,i-1))/dx+T2(1,i)*(P2(1,i)-P2(1,i-1))/(P2(1,i)*dx));
        T3(1,i)=-u2(1,i)*(T2(1,i)-T2(1,i-1))/dx-(gama-1)*T2(1,i)*((u2(1,i)-u2(1,i-1))/dx+u2(1,i)*(log(A(1,i))-log(A(1,i-1)))/dx);
        %计算导数平均值（7-60，7-61，7-62）
        P4(1,i)=0.5*(P1(1,i)+P3(1,i));
        u4(1,i)=0.5*(u1(1,i)+u3(1,i));
        T4(1,i)=0.5*(T1(1,i)+T3(1,i));
        %计算校正值（7-63，7-64，7-65）
        P(1,i)=P(1,i)+dt*P4(1,i);
        u(1,i)=u(1,i)+dt*u4(1,i);
        T(1,i)=T(1,i)+dt*T4(1,i);
        Ma(1,i)=u(1,i)/sqrt(T(1,i));
    end 
    %外流场计算（7-72）
    u(1,Nx)=2*u(1,Nx-1)-u(1,Nx-2);
    P(1,Nx)=2*P(1,Nx-1)-P(1,Nx-2);
    T(1,Nx)=2*T(1,Nx-1)-T(1,Nx-2);       
end
%显示结果
plot(x,P,'Linewidth',1.2,'Color','b'); hold on;
disp(P);
legend("P")


         
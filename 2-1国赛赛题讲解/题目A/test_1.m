clear
global m r g k F l
% m 铁锅质量
% r 铁锅半径
% g 重力加速度
% k 弹簧倔强系数
% F A端施加的向下的力的大小
% l 弹簧原长
m=5; r=0.3; g=9.8; k=500; F=50; l=0.8; T=20; % T:时长
[x,v]=fsolve(@fit,[0,0,-1,0,0,0]);
a=x(4);        %0<a<pi/2
b=x(5);        %0<b<2*pi,   n3 方向角
c=x(6);        % 确定n1方向  0<c<2*pi
O_0=x(1:3);
n3_0=[sin(a)*cos(b),sin(a)*sin(b),cos(a)];
tmp=[cos(a)*cos(b),cos(a)*sin(b),-sin(a)];   %tmp垂直与n3
tmp1=cross(n3_0,tmp);                         % tmp1,tmp2 锅面的两个垂直方向
n1_0=tmp*cos(c)+tmp1*sin(c);                  % n1方向
n2_0=cross(n3_0,n1_0);    
A=O_0+r*n1_0;
B=O_0-r/2*n1_0+sqrt(3)/2*r*n2_0;            %A,B,C,D坐标
C=O_0-r/2*n1_0-sqrt(3)/2*r*n2_0;
D=O_0-r*cos(2*pi/9)*n1_0+sin(2*pi/9)*r*n2_0;  
theta=acos(n3_0(3))/pi*180;


H1=5*m*r^2/12;
H2=H1;
H3=2*m*r^2/3;              %主轴方向的转动惯量（关于重心）

h=0.001;           %时间步长
n=floor(T/h)+1;
t=h*(0:(n-1))';
v=zeros(n,3);      %速度
n3=v;              %锅面法向
n1=v;              %n1方向
n2=v;              %n2方向
O=v;               %球心坐标
w=v;               %角速度
v(1,:)=[0,0,0];
n3(1,:)=n3_0;
n1(1,:)=n1_0;
n2(1,:)=n2_0;
O(1,:)=O_0;
w(1,:)=[0,0,0];
for i=1:n-1
    A=O(i,:)+r*n1(i,:);
    B=O(i,:)-r/2*n1(i,:)+sqrt(3)/2*r*n2(i,:);            %A,B,C,D坐标
    C=O(i,:)-r/2*n1(i,:)-sqrt(3)/2*r*n2(i,:);
    p=O(i,:)-r/2*n3(i,:);                                %重心p坐标

    
    l1=A-p;
    l2=B-p;
    l3=C-p;                               %弹簧拉力的力臂
    
    
    e1=-A/norm(A);
    e2=-B/norm(B);
    e3=-C/norm(C);                                  %弹簧A，B,C的单位方向（向上）
    
    F1=k*(norm(A)-l)*e1;
    F2=k*(norm(B)-l)*e2;
    F3=k*(norm(C)-l)*e3;                             %弹簧拉力
    M=cross(l1,F1)+cross(l2,F2)+cross(l3,F3);        %合力矩大小
    
    M1=dot(M,n1(i,:));
    M2=dot(M,n2(i,:));
    M3=dot(M,n3(i,:));                                    %力矩在主轴上的分量
    
    
    vc=v(i,:)+((F1+F2+F3)/m-g*[0,0,1])*h;
    Oc=O(i,:)+v(i,:)*h;
    
    w1c=w(i,1)+M1/H1*h;
    w2c=w(i,2)+M2/H2*h;
    w3c=w(i,3)+M3/H3*h;
    
    n1c=n1(i,:)+(-w(i,2)*n3(i,:)+w(i,3)*n2(i,:))*h;
    n2c=n2(i,:)+(-w(i,3)*n1(i,:)+w(i,1)*n3(i,:))*h;
    n3c=n3(i,:)+(-w(i,1)*n2(i,:)+w(i,2)*n1(i,:))*h;
    n1c=n1c/norm(n1c);
    n2c=n2c/norm(n2c);
    n3c=n3c/norm(n3c);
    
    A=Oc+r*n1c;
    B=Oc-r/2*n1c+sqrt(3)/2*r*n2c;            %A,B,C坐标
    C=Oc-r/2*n1c-sqrt(3)/2*r*n2c;
    p=Oc-r/2*n3c;                                %重心p坐标

    
    l1=A-p;
    l2=B-p;
    l3=C-p;                               %弹簧拉力的力臂
    
    
    e1=-A/norm(A);
    e2=-B/norm(B);
    e3=-C/norm(C);                                  %弹簧A，B,C的单位方向（向上）
    
    F1=k*(norm(A)-l)*e1;
    F2=k*(norm(B)-l)*e2;
    F3=k*(norm(C)-l)*e3;                             %弹簧拉力
    M=cross(l1,F1)+cross(l2,F2)+cross(l3,F3);        %合力矩大小
    
    M1=dot(M,n1c);
    M2=dot(M,n2c);
    M3=dot(M,n3c);                                    %力矩在主轴上的分量

    vp=v(i,:)+((F1+F2+F3)/m-g*[0,0,1])*h;
    Op=O(i,:)+vc*h;
    
    w1p=w(i,1)+M1/H1*h;
    w2p=w(i,2)+M2/H2*h;
    w3p=w(i,3)+M3/H3*h;
    
    n1p=n1(i,:)+(-w2c*n3c+w3c*n2c)*h;
    n2p=n2(i,:)+(-w3c*n1c+w1c*n3c)*h;
    n3p=n3(i,:)+(-w1c*n2c+w2c*n1c)*h;
    n1p=n1p/norm(n1p);
    n2p=n2p/norm(n2p);
    n3p=n3p/norm(n3p);
    
    O(i+1,:)=(Oc+Op)/2;
    v(i+1,:)=(vc+vp)/2;
    w(i+1,:)=[w1c+w1p,w2c+w2p,w3c+w3p]/2;
    
    n1(i+1,:)=(n1c+n1p)/2;
    n2(i+1,:)=(n2c+n2p)/2;
    n3(i+1,:)=(n3c+n3p)/2;
    n1(i+1,:)=n1(i+1,:)/norm(n1(i+1,:));
    n2(i+1,:)=n2(i+1,:)/norm(n2(i+1,:));
    n3(i+1,:)=n3(i+1,:)/norm(n3(i+1,:));
end

A=O+r*n1;
B=O-r/2*n1+sqrt(3)/2*r*n2;            %A,B,C坐标
C=O-r/2*n1-sqrt(3)/2*r*n2;
p=O-r/2*n3;                          %重心坐标

la=sqrt(sum(A.*A,2));
lb=sqrt(sum(B.*B,2));              %弹簧长度
lc=sqrt(sum(C.*C,2));

theta1=acos(abs(A(:,3))./la)*180/pi;
theta2=acos(abs(B(:,3))./lb)*180/pi;      %弹簧倾角
theta3=acos(abs(C(:,3))./lc)*180/pi;

theta=acos(n3(:,3))*180/pi;        %锅面倾角(度）

e=k*((la-l).^2+(lb-l).^2+(lc-l).^2)/2+m*g*p(:,3)+m*sum(v.*v,2)/2+(H1*w(:,1).^2+H2*w(:,2).^2+H3*w(:,3).^2)/2;
figure
ax1=subplot(4,2,1);
ax2=subplot(4,2,2);
ax3=subplot(4,2,3);
ax4=subplot(4,2,4);
ax5=subplot(4,2,5);
ax6=subplot(4,2,6);
ax7=subplot(4,2,7);
ax8=subplot(4,2,8);
plot(ax1,t,la-l);
title(ax1,'A端弹簧伸长量（米）');

plot(ax2,t,lb-l);
title(ax2,'B端弹簧伸长量（米）');

plot(ax3,t,lc-l)
title(ax3,'C端弹簧伸长量（米）');

plot(ax4,t,theta);
title(ax4,'锅面倾斜角（度）');

plot(ax5,t,theta1);
title(ax5,'A端弹簧与垂直方向夹角（度）');

plot(ax6,t,theta2);
title(ax6,'B端弹簧与垂直方向夹角（度）');

plot(ax7,t,theta3);
title(ax7,'C端弹簧与垂直方向夹角（度）');

plot(ax8,t, e);
title(ax8,'系统总能量（动能+势能）');




function y=fit(x)
global m r g k F l
% x(1:3) O点坐标
a=x(4);        %0<a<pi/2
b=x(5);        %0<b<2*pi,   n3 方向角
c=x(6);        % 确定n1方向  0<c<2*pi
O=x(1:3);
n3=[sin(a)*cos(b),sin(a)*sin(b),cos(a)];
tmp=[cos(a)*cos(b),cos(a)*sin(b),-sin(a)];   %tmp垂直与n3
tmp1=cross(n3,tmp);                         % tmp1,tmp2 锅面的两个垂直方向
n1=tmp*cos(c)+tmp1*sin(c);                  % n1方向
n2=cross(n3,n1);    
A=O+r*n1;
B=O-r/2*n1+sqrt(3)/2*r*n2;            %A,B,C,D坐标
C=O-r/2*n1-sqrt(3)/2*r*n2;
D=O-r*cos(2*pi/9)*n1+sin(2*pi/9)*r*n2;
p=O-r/2*n3;                          %重心p坐标


e1=-A/norm(A);
e2=-B/norm(B);
e3=-C/norm(C);  %弹簧A，B,C的单位方向（向上）

F1=k*(norm(A)-l)*e1;
F2=k*(norm(B)-l)*e2;
F3=k*(norm(C)-l)*e3;                             %弹簧A，B,C拉力
y(1:3)=F1+F2+F3-(F+m*g)*[0,0,1];                 % 合力大小

l1=A-p;
l2=B-p;
l3=C-p;                                             %弹簧拉力的力臂
l4=D-p;

y(4:6)=cross(l1,F1)+cross(l2,F2)+cross(l3,F3)-F*cross(l4,[0,0,1]);  % 合力矩大小
end






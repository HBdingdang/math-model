global dt
n=8;L=1.7; R=1.88; h=0.11; r=0.2; a=0.4; e=1;
F1=20; F2=80;
V1=zeros(50,1); i=1;
dt=0.001;
while 1
    tmp=(F1+F2)/2;  
    [~,speed,height]=drum(n,L,R,tmp);
    collision_height=height(length(height));
    [ ~,~,v2 ] = balldescend(a+h,collision_height );
    v=v2(length(v2));
    V1(i)=v;i=i+1;
    v1=collision(speed,v,e);
    if v1>min_ascend_speed(a+h,collision_height)
        F2=tmp;
    else
        F1=tmp;
    end
    if F2-F1<=1e-4
        break;
    end
end

F=(F1+F2)/2;
[time,speed,height]=drum(n,L,R,F);
collision_height=height(length(height));
t0=time(length(time));
figure
hold on
[ t,pos,v2 ] = balldescend(a+h,collision_height );
v=v2(length(v2));
t1=t(length(t));
plot(t,pos,'b');
plot(t1-t0+time,height,'r');
H=zeros(100,1);H(1)=a; i=2;
while 1
    if t1>30
        break;
    end
    tmp=collision(speed,v,e);
    [ t,pos ] = ballascend( collision_height,tmp);
    H(i)=pos(length(pos))-h;i=i+1;
    plot(t1+t,pos,'b');
    t1=t1+t(length(t));
    [t,pos,v2]=balldescend(pos(length(pos)),collision_height);
    v=v2(length(v2));
    plot(t1+t,pos,'b');
    t1=t1+t(length(t));
    plot(t1-t0+time,height,'r');
end

legend('������˶��켣','ͬ�Ĺ��������켣');
xlabel('ʱ�䣨�룩');
ylabel('�߶ȣ��ף�');
title(['������Ϊ',num2str(n),'������Ϊ',num2str(L),'�ס��˾������Ϊ',num2str(R),'�ס�������ĵĳ�ʼ����Ϊ',num2str(a),'�ס���ײ�ָ�ϵ��Ϊ',num2str(e),'������Ϊ',num2str(F),'ţ��ʱ��������ͬ�Ĺĵ��˶��켣']);


function [time,speed,height ] = drum( n,L,R,F )
global dt
% ͬ�Ĺĵ��˶����� 
% n: ��ҵĸ���
% L: ����
% R: ��Ҿ���ĵľ���
% F: ���ϵ�����
% time: ͬ�Ĺ��ٶ����ʱ����ʱ����
% speed: ͬ�Ĺĵ�����ٶ�
% height: ͬ�ĹĵĹ�������Ӧʱ�����ϵĸ߶�
h=0.11;
t=0:dt:10;
m=length(t);
y=zeros(m,2); y(1,:)=[0,0];
for i=1:m-1
    tmp=odedrum(y(i,:),n,L,R,F)*dt;
    yc=y(i,:)+tmp';
    tmp=odedrum(y(i,:),n,L,R,F)*dt;
    yp=y(i,:)+tmp';
    y(i+1,:)=(yp+yc)/2;
    if y(i+1,2)<y(i,2)
        time=t(1:i);speed=y(i,2);height=y(1:i,1)+h;
        break;
    end  
end
end

function dydt=odedrum(y,n,L,R,F)
% y(1) ���������꣬y(2) ���Ĵ�ֱ�����ٶ�
dydt=zeros(2,1);
M=3.6; r=0.2;C=1;rho=1.2258;g=9.8; 
H0=sqrt(L^2-(R-r)^2); % �Ĵ�ƽ��λ�õ���ˮƽλ�õľ���
S=pi*r^2;            % ͬ�Ĺ����
t1=n*F/M/L;  t2=C*rho*S/2/M;
dydt(1)=y(2);
dydt(2)=t1*(H0-y(1))-t2*y(2)^2-g;
end

function dydt=balldrop(y)
% y(1) ���������꣬y(2) ��ֱ�����ٶ�
dydt=zeros(2,1);
C=0.5;rho=1.2258; m=0.27;r=0.21; g=9.8;
s=pi*r^2;
t1=C*rho*s/2/m;
dydt(1)=y(2);
dydt(2)=-g+t1*y(2)^2;
end

function dydt=ballup(y)
% y(1) ���������꣬y(2) ��ֱ�����ٶ�
dydt=zeros(2,1);
C=0.5;rho=1.2258; m=0.27;r=0.21; g=9.8;
s=pi*r^2;
t1=C*rho*s/2/m;
dydt(1)=y(2);
dydt(2)=-g-t1*y(2)^2;
end


function [ time,pos,v ] = balldescend( a,b )
%  a: ����ʼ�߶�
%  b: �������ĸ߶�
%  time: ���a�䵽b����ʱ����
%  pos: ������Ӧʱ�����ϵĸ߶�
%  v:������Ӧʱ�����ϵ��ٶ�
global dt
t=0:dt:10;
n=length(t);
y=zeros(n,2);
y(1,:)=[a,0];
for i=1:n-1
    if y(i,1)<b
        time=t(1:i);
        pos=y(1:i,1);
        v=y(1:i,2);
        break;
    end
    tmp=balldrop(y(i,:))*dt;
    yc=y(i,:)+tmp';
    tmp=balldrop(yc)*dt;
    yp=y(i,:)+tmp';
    y(i+1,:)=(yc+yp)/2;
end
end

function v = collision( v1,v2 ,e)
% �ġ�����ײ����ٶ�
% v1: ����ײǰ��ֱ�����ٶ�  +�����ٶ�����
% v2: ����ײǰ��ֱ�����ٶ�  
% v: ����ײ��ֱ�����ٶ�  
% e: �ǵ�����ײ�ָ�ϵ�� e=1,������ײ
% M: ͬ�Ĺĵ�����
% m: �������
M=3.6;m=0.27;
if v1>v2
    v=((m-e*M)*v2+M*(1+e)*v1)/(m+M);
else
    v=v2;
end
end

function speed = min_ascend_speed(a, b )
global dt
%ascend�� ����������
%  a: ���������С����
%  b: ����ײ�߶�
% speed: ʹ��������ĸ߶ȳ���a�׵���С�ٶ�
t=0:-dt:-10;
n=length(t);
y=zeros(n,2);
for i=1:n-1
    if y(i,1)<=b-a
        speed=y(i,2);
        break;
    end
    tmp=-dt*ballup(y(i,:));
    yc=y(i,:)+tmp';
    tmp=-dt*ballup(yc);
    yp=y(i,:)+tmp';
    y(i+1,:)=(yc+yp)/2;
end
end

function [ time,height ] = ballascend( b,v0 )
%ascend�� ����������
%  b: ����ʼ�߶�
%  v0:����ʼ�ٶ�
%  time: ����������ߵ�����ʱ����
%  height: ������Ӧʱ�����ϵĸ߶�ֵ
global dt
t=0:dt:10;
n=length(t);
y=zeros(n,2);
y(1,:)=[b,v0];

for i=1:n-1
    if y(i,2)<0
        time=t(1:i);
        height=y(1:i,1);
        break;
    end
    tmp=ballup(y(i,:))*dt;
    yc=y(i,:)+tmp';
    tmp=ballup(yc)*dt;
    yp=y(i,:)+tmp';
    y(i+1,:)=(yp+yc)/2;
end
end

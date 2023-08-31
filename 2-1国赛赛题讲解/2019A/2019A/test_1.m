global m r g k F
% m ��Բ�α�Ƭ����
% r ��Ƭ�뾶
% g �������ٶ�
% k ���ɾ�ǿϵ��
% F B��ʩ�ӵ����µ����Ĵ�С
% ����ԭ��Ϊ����ԭ��λ��
m=5; r=0.3; g=9.8; k=500; F=50; T=600; % T:ʱ��
L=4*r/3/pi;    % ��Ƭ����λ��
y=(m*g+F)/k;      %�������ɵ����쳤
y1=0;         
y2=y/2;
H=m*(r^2/2+L^2+8*r*L/3/pi);  %ת������
dt=0.001;
n=floor(T/dt)+1;
t=dt*(0:(n-1))';
while 1
    tmp=(y1+y2)/2;
    s=fun(tmp,y-tmp);
    if s>=0
        y1=tmp;
    else 
        y2=tmp;
    end
    if y2-y1<1e-8
        y1=(y1+y2)/2;
        break;
    end
end
t1=(y-2*y1)/2/r;
theta=asin(t1);  %��ת�Ƕ�
y0=y/2+sqrt(1-t1^2)*L;  %����λ��

y0=[-y0,0,theta,0];
y=zeros(n,4);
y(1,:)=y0;
tic
for i=1:n-1
    tmp=odefcn(t,y(i,:))*dt;
    yc=y(i,:)+tmp';
    tmp=odefcn(t,yc)*dt;
    yp=y(i,:)+tmp';
    y(i+1,:)=(yc+yp)/2;    
end
toc
theta=y(:,3);
w=y(:,4);
v=y(:,2);
x3=y(:,1);
x1=-x3-L*cos(theta)-r*sin(theta);
x2=-x3-L*cos(theta)+r*sin(theta);
e=1/2*k*(x1.^2+x2.^2)+1/2*H*w.^2+1/2*m*v.^2+m*g*x3;

figure
ax1=subplot(4,2,1);
ax2=subplot(4,2,2);
ax3=subplot(4,2,3);
ax4=subplot(4,2,4);
ax5=subplot(4,2,5);
ax6=subplot(4,2,6);
ax7=subplot(4,2,7);
plot(ax1,t,x1);
title(ax1,'A�˵����쳤�����ף�(�Ľ�ŷ����)');
plot(ax2,t,x2);
title(ax2,'B�˵����쳤�����ף�(�Ľ�ŷ����)');
plot(ax3,t,x3);
title(ax3,'��Բ�����������꣨�ף�');
plot(ax4,t,v);
title(ax4,'��Բ�̴�ֱ�����ٶȣ���/�룩');
plot(ax5,t,theta);
title(ax5,'��Բ����б�ǣ����ȣ�');
plot(ax6,t,w);
title(ax6,'��Բ�������ĵ���ת���ٶȣ�����/�룩');
plot(ax7,t, e);
title(ax7,'ϵͳ������������+���ܣ�');


tic
[~,y]=ode45(@(t,y)odefcn(t,y),t,y0');
toc
theta=y(:,3);
w=y(:,4);
v=y(:,2);
x3=y(:,1);
x1=-x3-L*cos(theta)-r*sin(theta);
x2=-x3-L*cos(theta)+r*sin(theta);
e=1/2*k*(x1.^2+x2.^2)+1/2*H*w.^2+1/2*m*v.^2+m*g*x3;

figure
ax1=subplot(4,2,1);
ax2=subplot(4,2,2);
ax3=subplot(4,2,3);
ax4=subplot(4,2,4);
ax5=subplot(4,2,5);
ax6=subplot(4,2,6);
ax7=subplot(4,2,7);
plot(ax1,t,x1);
title(ax1,'A�˵����쳤�����ף�(ode45)');
plot(ax2,t,x2);
title(ax2,'B�˵����쳤�����ף�(ode45)');
plot(ax3,t,x3);
title(ax3,'��Բ�����������꣨�ף�');
plot(ax4,t,v);
title(ax4,'��Բ�̴�ֱ�����ٶȣ���/�룩');
plot(ax5,t,theta);
title(ax5,'��Բ����б�ǣ����ȣ�');
plot(ax6,t,w);
title(ax6,'��Բ�������ĵ���ת���ٶȣ�����/�룩');
plot(ax7,t, e);
title(ax7,'ϵͳ������������+���ܣ�');




function dydt = odefcn(t,y)
global  r m g k
L=4*r/3/pi;  
H=m*(r^2/2+L^2+8*r*L/3/pi);
y1=-y(1)-L*cos(y(3))-r*sin(y(3));
y2=-y(1)-L*cos(y(3))+r*sin(y(3));
x1=L*sin(y(3))-r*cos(y(3));
x2=L*sin(y(3))+r*cos(y(3));
dydt=zeros(4,1);
dydt(1)=y(2);
dydt(2)=k*(y1+y2)/m-g;
dydt(3)=y(4);
dydt(4)=-k*(x1*y1+x2*y2)/H;
end


function v=fun(y1,y2)
%y1    A�˵����쳤
%y2    B�˵����쳤
global  r  k F
L=4*r/3/pi;    
stmp=(y2-y1)/2/r;  %ˮƽ��ǵ�����ֵ
ctmp=sqrt(1-stmp^2);
x0=L*stmp;
x2=x0+r*ctmp;
x1=x0-r*ctmp;
v=x2*(k*y2-F)+x1*k*y1;    %�������ĵ����ش�С
end
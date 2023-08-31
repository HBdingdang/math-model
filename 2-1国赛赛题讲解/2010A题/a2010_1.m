clc,clear;  %����λ��ģ��
global a b L;
a=1.78/2; b=1.2/2; L=2.45; alpha1=0;alpha2=4.1/180*pi;

filename="����A����1��ʵ��ɼ����ݱ�.xls";
A=xlsread(filename,1,'C:C')+262;  % ʵ������ ��λ�� ��δ��
B=xlsread(filename,1,'D:D');      % ����߶� ��λ mm
A1=xlsread(filename,3,'C:C')+215;  % ʵ������ ��λ�� ������λ��4.1��
B1=xlsread(filename,3,'D:D');

v0=zeros(size(B));
for i=1:length(B)
    v0(i)=f3(B(i)/1000,alpha1); %������������
end    
v1=zeros(size(B1));
for i=1:length(B1)
    v1(i)=f3(B1(i)/1000,alpha2); %������������
end

tmp1=zeros(1200,1); tmp2=tmp1;
for i=1:1200
    tmp1(i)=f3(i/1000,alpha1);
    tmp2(i)=f3(i/1000,alpha2);
end

figure
ax1=subplot(2,1,1); ax2=subplot(2,1,2); 
plot(ax1,1:1200,tmp1,1:1200,real(tmp2));
title(ax1,'��λ��(����������λ��4.1��ʱ����������');
plot(ax2,B,A, B1,A1);
title(ax2,'��λ�䣨����������λ��4.1��ʱ��ʵ������');
c0=v0-A;c1=v1-A1;
figure
ax1=subplot(2,2,1);
plot(ax1,B,v0,B,A,'--');   %���Կ�������������ʵ��������ƫ��
title(ax1,'��������(������ʵ����������λ��')
ax2=subplot(2,2,3);
plot(ax2,B,c0);     %����������ʵ��������ƫ��
title(ax2,'����������ʵ�������Ĳ�ֵ����λ��');

ax3=subplot(2,2,2);
plot(ax3,B1,v1,B1,A1,'--');   %���Կ�������������ʵ��������ƫ��
title(ax3,'��������(��)��ʵ������������λ��4.1��')
ax4=subplot(2,2,4);
plot(ax4,B1,c1);     %����������ʵ��������ƫ��
title(ax4,'����������ʵ�������Ĳ�ֵ������λ��4.1��');


c2=zeros(size(B1));
for i=1:length(B1)
    c2(i,1)=f4(B1(i));
end
tmp=c1-c2;
c=sum(tmp)/length(tmp);

%�ɳ�������֤���
A3=xlsread(filename,2,'C:C'); % �ۼӳ����� ��λ�� ��δ��
B3=xlsread(filename,2,'D:D');      % ����߶� ��λ mm

t1=f3(B3(1)/1000,alpha1)-f4(B3(1))+A3(1);

v3=zeros(size(B3));
for i=1:length(B3)
    v3(i)=t1-f3(B3(i)/1000,alpha1)+f4(B3(i));
end


A4=xlsread(filename,4,'C:C');  % �ۼӳ����� ��λ�� ����λ��4.1��
B4=xlsread(filename,4,'D:D');      % ����߶� ��λ mm

%������������֤��ʹ��λ����
t2=f3(B4(1)/1000,alpha2)-f5(B4(1))+A4(1);

v4=zeros(size(B4));
for i=1:length(B4)
    v4(i)=t2-f3(B4(i)/1000,alpha2)+f5(B4(i));
end

tmp1=f3(B4(1)/1000,alpha2)-f5(B4(1));

v5=zeros(length(B4)-1,1);
for i=1:length(B4)-1
    tmp2=f3(B4(i+1)/1000,alpha2)-f5(B4(i+1));
    v5(i)=tmp1-tmp2;
    tmp1=tmp2;
end

%������������֤��ʹ���ޱ�λ��
t3=f3(B4(1)/1000,alpha2)-f4(B4(1))-c+A4(1);

v6=zeros(size(B4));
for i=1:length(B4)
    v6(i)=t3-f3(B4(i)/1000,alpha2)+f4(B4(i))+c;
end

tmp1=f3(B4(1)/1000,alpha2)-f4(B4(1));
v7=zeros(length(B4)-1,1);
for i=1:length(B4)-1
    tmp2=f3(B4(i+1)/1000,alpha2)-f4(B4(i+1));
    v7(i)=tmp1-tmp2;
    tmp1=tmp2;
end

tmp1=v3-A3;tmp2=v4-A4; tmp3=v6-A4;%�����ۼӳ�������ʵ���ۼӳ������Ĳ�ֵ

figure
ax1=subplot(5,1,1);ax2=subplot(5,1,2);ax3=subplot(5,1,3);ax4=subplot(5,1,4);ax5=subplot(5,1,5);
plot(ax1,B3,tmp1);   
title(ax1,'�����ۼӳ�������ʵ���ۼӳ�������ֵ����λ��');

plot(ax2,B4,tmp2);
title(ax2,'�����ۼӳ�������ʵ���ۼӳ�������ֵ(ʹ�ñ�λ������)������н�4.1��');

plot(ax3,B4,tmp3);
title(ax3,'�����ۼӳ�������ʵ���ۼӳ�������ֵ��ʹ���ޱ�λ�����㣩������н�4.1��');

plot(ax4,1:length(B4)-1,v5);
title(ax4,'��λ�������������ÿ�γ�������50��������н�4.1��');

plot(ax5,1:length(B4)-1,v7);
title(ax5,'����λ�������������ÿ�γ�������50��������н�4.1��');

%�����궨4.1�������
youliang=zeros(121,1);
for i=0:120
    youliang(i+1)=f3(i/100,alpha2)-f4(10*i)-c;
end

figure
plot(0:10:1200,youliang,B1,A1);
title('�궨��������ʵ��������4.1�㣩');

function s=f(h) %��h����Բ�������
global a b;
  s=a*b.*acos((b-h)/b)-(b-h).*a.*sqrt(1-(h-b).^2/b^2);
end


function l=f1(h,alpha) %��������Ϊhʱ����ֱ��h���������߶Ⱥ��� x(h)
   global L;
   h1=2.05*tan(alpha);
   if alpha==0
       l=L;
   elseif h<=h1
       l=0.4+cot(alpha)*h;
   else 
       l=L;
   end
end


function height=f2(x,h,alpha) %��������λhʱ����x������Բ����߶�   h(l)
    global b;
    if alpha==0
        height=h*ones(size(x));
    else
        h2=2*b-0.4*tan(alpha);
        tmp=0.4-(2*b-h)*cot(alpha);
        if h<=h2
            height=h+(0.4-x)*tan(alpha);
        else
            height=2*(x<=tmp).*b+(x>tmp).*(h+(0.4-x).*tan(alpha));
        end
    end
end

function v=f3(h,alpha)  %������hʱ����������
    v=integral(@(x) f(f2(x,h,alpha)),0,f1(h,alpha))*1000;
end


%��cftool��Ϲ��ߵ��ݻ�
function p=f4(h)  %����ƫ��� h ��λ(mm), ��λ��
    p=-8.403e-08*h.^3+0.0001506*h.^2+0.05822*h-1.711;
end

function p=f5(h)  %����ƫ��� ��λ��alpha=4.1
   p=2.881e-07*h.^2-0.001026*h.^2+1.023*h-222.1;
end
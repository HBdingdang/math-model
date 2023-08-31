clc,clear;  %����λ��ģ��
global R r L ;
r=1.5; R=13/8; L=8;
filename="����A����2��ʵ�ʲɼ����ݱ�.xls";
A=xlsread(filename,1,'E:E'); % ��λ�߶�mm
B=xlsread(filename,1,'F:F'); % ��ʾ��λ����L
n=size(A); m=302; C=zeros(n);
for i=1:n
   C(i)=volume1(A(i)/1000,0,0)+volume2(A(i)/1000,0,0)+volume3(A(i)/1000,0,0); %��������
end
D=B-C; [a,I]=sort(A);
d=D(I);

figure
ax1=subplot(5,1,1);ax2=subplot(5,1,2);ax3=subplot(5,1,3); ax4=subplot(5,1,4); ax5=subplot(5,1,5);
plot(ax1, A,B,A,C,'--');
title(ax1,'�ޱ�λ�����µ�ʵ������ֵ�����ۼ���ֵ');

plot(ax2,A,D);
title(ax2,'�ޱ�λ�����µ�ʵ������ֵ�����ۼ���ֵ֮��');

plot(ax3,A(1:m),D(1:m));
title(ax3,'�ޱ�λ�����µ�ʵ�����������ۼ���ֵ�Ĳ�ֵ��1-302��');

plot(ax4,A(m+1:n),D(m+1:n));
title(ax4,'�ޱ�λ�����µ�ʵ�����������ۼ���ֵ�Ĳ�ֵ��303-603��');

plot(ax5,a,d);
title(ax5,'�ޱ�λ����������ĵ�ʵ�����������ۼ���ֵ�Ĳ�ֵ');


%��λ�궨����alpha=2.1��,beta=4.4�㣩
i=10:10:2*r*100; c=zeros(size(i));
for j=1:length(i)
    c(j)=volume(10*i(j),2.1*pi/180,4.4*pi/180);
    c1(j)=volume(10*i(j),0,0);
end
    
figure
ax1=subplot(1,1,1);
plot(ax1,a,B(I),'*',i*10,c1,'b--o',i*10,c,i*10,c-c1,'--');
title(ax1,'ԭ��λ����(*)����������(-o-)��ʵ������(��)����ֵ�Ա�ͼ���������2.1��,�������4.4��');
    

E=xlsread(filename,1,'D:D'); % ÿ�γ���������

%��ʵ��������alpha,beta ȡ��һ�γ��������� ���� alpha=2.1, beta=4.4, M=247.988
k=1:302;  
alpha_0=2.1;beta_0=4.4;M=0; 
tmp=volume(A(k(1)),alpha_0*pi/180,beta_0*pi/180);
for s=1:length(k)-1
    tmp1=volume(A(k(s+1)),alpha_0*pi/180,beta_0*pi/180);
    M=M+(tmp-tmp1-E(k(s+1)))^2;
    tmp=tmp1;
end

for i=1:101
    for j=1:101
        M1=0; tmp=volume(A(k(1)),(i-1)*pi/1800,(j-1)*pi/1800);
        for s=1:length(k)-1
            tmp1=volume(A(k(s+1)),(i-1)*pi/1800,(j-1)*pi/1800);
            M1=M1+(tmp-tmp1-E(k(s+1)))^2;
            tmp=tmp1;
            if M1>=M
                break;
            end
        end
        if M1<M
            alpha_0=(i-1)/10;beta_0=(j-1)/10;M=M1;
        end
    end
end

%�ɵڶ�������304--603���� ���������֤ ���Ϊ293.2904
k=304:603;  
alpha_0=2.1;beta_0=4.4;M=0; 
tmp=volume(A(k(1)),alpha_0*pi/180,beta_0*pi/180);
for s=1:size(k)-1
    tmp1=volume(A(k(s+1)),alpha_0*pi/180,beta_0*pi/180);
    M=M+(tmp-tmp1-E(k(s+1)))^2;
    tmp=tmp1;
end
    
           
        
% ���岿��

function s=section_area(r1,h) 
  s=(h-r1).*sqrt(r1.^2-(r1-h).^2)+r1.^2.*acos(1-h./r1);
end

function l=length1(h,alpha,beta) %���߶�Ϊhʱ����ֱ��h��������泤�Ⱥ��� x(h)
   global L r;
   h1=6*tan(alpha)*sec(beta)-r*(sec(beta)-1);
   if alpha==0
       l=L;
   elseif h<h1
       l=(h-r)*cos(beta)*cot(alpha)+r*cot(alpha)+2;
   else
       l=L;
   end
end


function height=height1(x,h,alpha,beta) %��������λhʱ����x������Բ����߶�  h(l)
    global r;
    if alpha==0
        height=((h-r)*cos(beta)+r)*ones(size(x));
    else
        tmp=2-(r-(h-r)*cos(beta))*cot(alpha);
        h2=r*(sec(beta)+1)-2*tan(alpha)*sec(beta);
        if h<=h2
            height=(h-r)*cos(beta)+r+(2-x)*tan(alpha);
        else
            height=(x<=tmp).*2*r+(x>tmp).*((h-r)*cos(beta)+r+(2-x).*tan(alpha));
        end
    end
end

function v=volume1(h,alpha,beta)  %������hʱ����������    ���岿�����
    global r;
    v=integral(@(x) section_area(r,height1(x,h,alpha,beta)),0,length1(h,alpha,beta))*1000;
end


%����ڲ���
function l=length2(h,alpha,beta)       %x1(h)
    global R r;
    tmp=(h-r)*cos(beta)+2*tan(alpha);
    l=(R-1+tmp*tan(alpha)-sqrt((1-R)^2-2*(1-R)*tmp*tan(alpha)-sec(alpha)^2*(1-2*R)-tmp^2))/sec(alpha)^2;
end

function height=height2(x,h,alpha,beta)  %h(x)
    global R r;
    h3=r-3*tan(alpha)*sec(beta);
    tmp=length2(h,alpha,beta);
    if h<h3
        height=(h-r)*cos(beta)-(x-2)*tan(alpha)+sqrt((2*R-x-1).*(x+1));
    else
        height=(x>=tmp).*(sqrt((2*R-x-1).*(x+1))+(h-r)*cos(beta)-(x-2)*tan(alpha))+(x<tmp).*(2*sqrt((2*R-x-1).*(x+1)));
    end
end

function v=volume2(h,alpha,beta)  %����ڲ������
    global R r;
    h3=r-3*tan(alpha)*sec(beta);
    h2=r*(sec(beta)+1)-2*tan(alpha)*sec(beta);
    if 2*tan(alpha)>=r*(1-cos(beta))&& h>=h2
        v=pi*integral(@(x) (2*R-x-1).*(x+1),-1,0)*1000;
    elseif h<h3
        v=integral(@(x) section_area(sqrt((2*R-x-1).*(x+1)),height2(x,h,alpha,beta)),length2(h,alpha,beta),0)*1000;
    else
        v=integral(@(x) section_area(sqrt((2*R-x-1).*(x+1)),height2(x,h,alpha,beta)),-1,0)*1000;
    end
end


%����ڲ��ֵ����
function l=length3(h,alpha,beta)
    global R r L;
    h1=6*tan(alpha)*sec(beta)-r*(sec(beta)-1);
    if h1>=0 && h<=h1
        l=L;
    else
        tmp=(h-r)*cos(beta)+2*tan(alpha);
        l=(9-R+tmp*tan(alpha)+sqrt((9-R)^2+2*(9-R)*tmp*tan(alpha)-9*sec(alpha)^2*(9-2*R)-tmp^2))/sec(alpha)^2;
    end
end

function height=height3(x,h,alpha,beta)   %x1(h)
    global L R r;
    h4=r+7*tan(alpha)*sec(beta);
    tmp=length3(h,alpha,beta);
    if h<=h4&&tmp>L
        height=(h-r)*cos(beta)-(x-2)*tan(alpha)+sqrt((2*R+x-9).*(9-x));
    elseif h<=h4&& tmp==L
        height=0;
    else
        height=(sqrt((2*R+x-9).*(9-x))+(h-r)*cos(beta)-(x-2)*tan(alpha)).*(x<=tmp)+(x>tmp).*(2*sqrt((2*R+x-9).*(9-x)));
    end
end

function v=volume3(h,alpha,beta)
    global R r;
    h1=6*tan(alpha)*sec(beta)-r*(sec(beta)-1);
    h4=r+7*tan(alpha)*sec(beta);
    if 6*tan(alpha)>=r*(1-cos(beta)) && h<=h1
        v=0;
    elseif h<=h4
        v=integral(@(x) section_area(sqrt((2*R+x-9).*(9-x)),height3(x,h,alpha,beta)),8,length3(h,alpha,beta))*1000;
    else
        v=integral(@(x) section_area(sqrt((2*R+x-9).*(9-x)),height3(x,h,alpha,beta)),8,9)*1000;
    end
end

function v=volume(h,alpha,beta) 
   v=volume1(h/1000,alpha,beta)+volume2(h/1000,alpha,beta)+volume3(h/1000,alpha,beta);
end


load('data.mat','l','L','P','E','e','index','side_ind_square');
global k R err1 err2  err3 p 
R=300.4;       %����뾶
f=0.466*R;     %����
s=-0.3361;
F=f-s;
alpha=36.795*pi/180;
beta=78.169*pi/180;
p=[cos(alpha)*cos(beta),cos(beta)*sin(alpha),sin(beta)];
H1=150^2/4/F*1.01^2;                    %�ھ���1%��������߶�
H0=4*(sqrt(F)*(H1+F)^(3/2)-F^2)/3/R;    %�����α������߶�
k=find(E(:,1:3)*p'<=H0-R);              %�α������ڵĽڵ�����
%l                    ����������ƽ��
%L                    �ٶ���ԭ��  
%e                    �ٶ�������
%P                    �ٶ����¶�����
%E                    �����ڵ�ԭ����
%side_ind_square      ������������������ƽ��
%index                �����йصľ���
%err1                 max(0,0.07%-�߳����������)  err1=0, ��ʾ���㲻��ʽ   
%err2                 ���������������
%err3                 �ٶ����������������  err3=0, |si|<0.6

err1=1; err2=1;err3=1;
load('x0.mat');
c=200000;    
g1=fun(x0);
g2=ineq(x0);
g3=eqf(x0);
g0=g1+c*(g2+g3);
r=2;                             %�Ŵ�ϵ��
eps=1e-6;                        %����
range=10;tol=1;
while tol>eps
     x1=x0;
     while 1
        y0=dfun(x0)+c*(dineq(x0)+deqf(x0));
        tmp=norm(y0);
        if tmp<eps
           break;
        end
        h=-y0/tmp;
        fhandle=@(t) fun(x0+t*h)+c*(ineq(x0+t*h)+eqf(x0+t*h));
        range=max(tmp,range);
        [t,g0]=funm(fhandle,range);
        range=t*50;
        g1=fun(x0);
        g2=ineq(x0);
        g3=eqf(x0);
        x0=x0+t*h;
    end
    if c>1e8
        break;
    else
        c=r*c;
    end
    tol=max([err1,err2,err3]); 
end



function y=dineq(x)                    %x(n,4) �����ڵ����꼰�ٶ���������
%y ����ʽԼ����x���ݶ�ֵ
global side_ind_square index           
r=0.0007;                               % �߳��������������
E=x(:,1:3);                             % �����ڵ�����
s=x(:,4);                               %�ٶ����쳤��
y=zeros(size(x));
tmp=E(side_ind_square(:,1),:)-E(side_ind_square(:,2),:);
tmp1=sum(tmp.*tmp,2);
tmp2=tmp1-side_ind_square(:,3)*(1+r)^2;
tmp3=side_ind_square(:,3)*(1-r)^2-tmp1;
tmp2=max(0,tmp2);
tmp3=max(0,tmp3);
tmp4=4*(tmp2-tmp3).*tmp;
tmp5=max(0,s.^2-0.36);
y(:,1:3)=index*tmp4;
y(:,4)=4*s.*tmp5;
end

function g=ineq(x)                    %x(n,4) �����ڵ����꼰�ٶ���������
%g ����ʽԼ����xֵ
global side_ind_square err1 err3 
r=0.0007;                               % �߳��������������
E=x(:,1:3);                             % �����ڵ�����
s=x(:,4);                               %�ٶ����쳤��
tmp=E(side_ind_square(:,1),:)-E(side_ind_square(:,2),:);
tmp1=sum(tmp.*tmp,2);
tmp2=tmp1-side_ind_square(:,3)*(1+r)^2;
tmp3=side_ind_square(:,3)*(1-r)^2-tmp1;
tmp2=max(0,tmp2);
tmp3=max(0,tmp3);
tmp4=max(0,s.^2-0.36);
g=sum(tmp2.^2)+sum(tmp3.^2)+sum(tmp4.^2);
err3=max(tmp4/0.72);
err1=max(max(tmp2,tmp3)./side_ind_square(:,3));
end

function g=eqf(x)       %x(n,4) �����ڵ����꼰�ٶ���������
%g ��ʽԼ����xֵ
global L l P e err2                       
E=x(:,1:3);                       % �����ڵ�����
s=x(:,4);                         %�ٶ����쳤��
tmp=P-E+e(:,1:3).*(L+s);
tmp1=sum(tmp.*tmp,2)-l;
t1=tmp1.^2;
g=sum(t1);
err2=max(abs(tmp1)./l)/2;
end




function y=deqf(x)       %x(n,4) �����ڵ����꼰�ٶ���������
%y ��ʽԼ����x���ݶ�ֵ
global L l P e                       
E=x(:,1:3);                       % �����ڵ�����
s=x(:,4);                         %�ٶ����쳤��
y=zeros(size(x));
tmp=P-E+e(:,1:3).*(L+s);
tmp1=sum(tmp.*tmp,2)-l;
tmp2=sum((P-E).*e,2);
y(:,1:3)=4*tmp1.*(E-P-(L+s).*e);
y(:,4)=4*tmp1.*(L+s+tmp2);
end

function g=fun(x)                  %x(n,4) �����ڵ����꼰�ٶ���������
%g �Ż�������xֵ
global k R p
f=0.466*R;                       
s=-0.3361;
F=f-s;                    %����
E=x(k,1:3);  
tmp=sum(E.^2,2);
tmp1=E(:,1:3)*p';
tmp2=tmp1.^2;
tmp3=tmp-tmp2-4*F*(tmp1-s+R);
g=tmp3'*tmp3;
end




function y=dfun(x)                  %x(n,4) �����ڵ����꼰�ٶ���������
%y �Ż�������x���ݶ�ֵ
global k R p
f=0.466*R;                       
s=-0.3361;
F=f-s;                    %����
y=zeros(size(x));
E=x(k,1:3);  
tmp=sum(E.^2,2);
tmp1=E(:,1:3)*p';
tmp2=tmp1.^2;
tmp3=tmp-tmp2-4*F*(tmp1-s+R);
y(k,1:3)=4*tmp3.*(E-tmp1*p-2*F*ones(length(k),1)*p);
end

function [x,y]=funm(f,x1)
 tmp1=f(0); tmp2=f(x1);
 if tmp1<=tmp2
     a=0;   fmin=tmp1; b=x1;
 else
     a=x1;b=0; fmin=tmp2;
 end
 while 1
    x=(a+b)/2;
    tmp=f(x);
    if tmp<fmin
        b=a; a=x; fmin=tmp;
    else
        b=x;
    end
    if abs(a-b)<1e-12
        x=a; y=fmin;
        break;
    end
 end      
end
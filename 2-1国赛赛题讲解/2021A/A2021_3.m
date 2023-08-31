clear
global R p
R=300.4;       %����뾶
f=0.466*R;     %����
s=-0.3361;
F=f-s;
alpha=36.795*pi/180;
beta=78.169*pi/180;
p=[cos(alpha)*cos(beta),cos(beta)*sin(alpha),sin(beta)];
H1=150^2/4/F;

load('data.mat','ind_area');
load('x0.mat');
E=x0(:,1:3);                           %�����ڵ�����
k=find(E*p'<=H1-R+s);                  %�����α�������ڵ�
S1=zeros(4300,2);
for i=1:4300
    A=ind_area(i,1);B=ind_area(i,2);C=ind_area(i,3);
    tag=sum(sum(k==[A,B,C]));           %�жϾ�������һ�����ڹ���������
    if tag
      A=E(A,:); B=E(B,:);C=E(C,:); u=ind_area(i,4);
      S1(i,:)=area1(A,B,C,u);
    end
end
rate1=sum(S1(:,1))/sum(S1(:,2));                %���ʷּ����300�׿ھ��Ļ�׼���淴���


load('data.mat','E');
H0=R-sqrt(R^2-150^2);                      %�ھ�300������߶�
k=find(E(:,3)<=H0-R);                      %���������ڵĽڵ�����
S2=zeros(4300,2);
for i=1:4300
    A=ind_area(i,1);B=ind_area(i,2);C=ind_area(i,3);
    tag=sum(sum(k==[A,B,C]));
    if tag
      A=E(A,:); B=E(B,:);C=E(C,:); u=ind_area(i,4);
      S2(i,:)=area2(A,B,C,u);
    end
end
rate=sum(S2(:,1))/sum(S2(:,2));                %���ʷּ����300�׿ھ��Ļ�׼���淴���


fun=@(x) -2*x.^2+(2.*x+1)*R^2-4*0.534*R*x.*sqrt(R^2-x.^2);          %ȷ����Ч��������ھ�
fun1=@(x) 2*x.^2+(2.*x-1)*R^2-4*0.534*R*x.*sqrt(R^2-x.^2);      
a=7; b=8;

while 1
    tmp=(a+b)/2;
    tmp1=fun(tmp);
    if tmp1>0
        a=tmp;
    else 
        b=tmp;
    end
    if b-a<1e-12
        r1=(a+b)/2;
        break;
    end
end

a=90; b=120;
while 1
    tmp=(a+b)/2;
    tmp1=fun(tmp);
    if tmp1>0
        b=tmp;
    else 
        a=tmp;
    end
    if b-a<1e-12
        r2=(a+b)/2;
        break;
    end
end

a=90; b=120;
while 1
    tmp=(a+b)/2;
    tmp1=fun1(tmp);
    if tmp1>0
        b=tmp;
    else 
        a=tmp;
    end
    if b-a<1e-12
        r3=(a+b)/2;
        break;
    end
end
rate2=(sqrt(R^2-0.25)-sqrt(R^2-r1^2)+sqrt(R^2-r2^2)-sqrt(R^2-r3^2))/(R-sqrt(R^2-150^2));    %���ڿھ�300�Ļ�׼�������Դ�ֽ��ܱ�

E=x0(:,1:3);                           %�����ڵ�����
k=find(E*p'<=H1-R+s);                    %�������ڵĽڵ�����
k1=zeros(4300,1);
for i=1:length(k)
    k1=k1+(ind_area(:,1)==k(i))+(ind_area(:,2)==k(i))+(ind_area(:,3)==k(i));
end
k1=find(k1~=0);                        %�������ڵľ�������(������һ�������ڹ�������)


n=1000;                             %���������
theta=rand(n,1)*2*pi;               %(0,2*pi)֮����ȷֲ����ݵ�
r=150*sqrt(rand(n,1));                  %(0,1)֮��r^2���ȷֲ����ݵ�    ������Ȳ���
n1=[cos(alpha)*sin(beta),sin(beta)*sin(alpha),-cos(beta)];
n2=cross(p,n1);
Qx=r.*cos(theta)*n1(1)+r.*sin(theta)*n2(1);
Qy=r.*cos(theta)*n1(2)+r.*sin(theta)*n2(2);
Qz=r.*cos(theta)*n1(3)+r.*sin(theta)*n2(3);        %�����������
Q=[Qx';Qy';Qz']';
counter=0;
A=E(ind_area(k1,1),:); B=E(ind_area(k1,2),:); C=E(ind_area(k1,3),:); 
a=B-A;b=C-A;
tmp=[a(:,2).*b(:,3)-b(:,2).*a(:,3),b(:,1).*a(:,3)-a(:,1).*b(:,3),a(:,1).*b(:,2)-a(:,2).*b(:,1)];                %������ƽ�淨��
c=sum(A.*tmp,2);
d=sum(tmp.*p,2);
for i=1:n
    t0=(c-sum(tmp.*Q(i,:),2))./d;
    D=Q(i,:)+t0.*p;
    e=D-A;
    tmp1=sum(a.*a,2);tmp2=sum(a.*b,2);tmp3=sum(b.*b,2);
    tmp4=sum(e.*a,2);tmp5=sum(e.*b,2);
    s=(tmp4.*tmp3-tmp5.*tmp2)./(tmp1.*tmp3-tmp2.*tmp2);
    t=(tmp1.*tmp5-tmp2.*tmp4)./(tmp1.*tmp3-tmp2.*tmp2);
    j=find((s>=0).*(t>=0).*(s+t<=1));
    if ~isempty(j)
        o=yuanxin(A(j,:),B(j,:),C(j,:));                         %��������
        tmp1=norm(D(j,:)-o);
        e=o+R/tmp1*(D(j,:)-o);                                %QD�뾵��Ľ���e
        q=2*(o-e)*p'*(o-e)/R^2-p;                               %������߷���
        tmp3=q*p';
        if(tmp3>0)
            t=(-0.534*R-e*p')/tmp3;
            tmp4=e+t*q+0.534*R*p;
            counter=counter+(norm(tmp4)<=0.5);
        end  
    end
end
rate3=counter/n;


%���Ȳ���
x=-150:3:150;y=x;
[x,y]=meshgrid(x,y);
r=sqrt(x.^2+y.^2);
k=find(r<=150);
r=r(k); x=x(k); y=y(k);
theta=(x>=0).*asin(y./(r+(r==0)))+(x<0).*(pi-asin(y./(r+(r==0))));
n1=[cos(alpha)*sin(beta),sin(beta)*sin(alpha),-cos(beta)];
n2=cross(p,n1);
Qx=r.*cos(theta)*n1(1)+r.*sin(theta)*n2(1);
Qy=r.*cos(theta)*n1(2)+r.*sin(theta)*n2(2);
Qz=r.*cos(theta)*n1(3)+r.*sin(theta)*n2(3);        %�����������
Q=[Qx(:)';Qy(:)';Qz(:)']';
counter=0;
A=E(ind_area(k1,1),:); B=E(ind_area(k1,2),:); C=E(ind_area(k1,3),:); 
a=B-A;b=C-A;
tmp=[a(:,2).*b(:,3)-b(:,2).*a(:,3),b(:,1).*a(:,3)-a(:,1).*b(:,3),a(:,1).*b(:,2)-a(:,2).*b(:,1)];                %������ƽ�淨��
c=sum(A.*tmp,2);
d=sum(tmp.*p,2);
for i=1:length(r)
    t0=(c-sum(tmp.*Q(i,:),2))./d;
    D=Q(i,:)+t0.*p;
    e=D-A;
    tmp1=sum(a.*a,2);tmp2=sum(a.*b,2);tmp3=sum(b.*b,2);
    tmp4=sum(e.*a,2);tmp5=sum(e.*b,2);
    s=(tmp4.*tmp3-tmp5.*tmp2)./(tmp1.*tmp3-tmp2.*tmp2);
    t=(tmp1.*tmp5-tmp2.*tmp4)./(tmp1.*tmp3-tmp2.*tmp2);
    j=find((s>=0).*(t>=0).*(s+t<=1));
    if ~isempty(j)
        o=yuanxin(A(j,:),B(j,:),C(j,:));                         %��������
        tmp1=norm(D(j,:)-o);
        e=o+R/tmp1*(D(j,:)-o);                                %QD�뾵��Ľ���e
        q=2*(o-e)*p'*(o-e)/R^2-p;                               %������߷���
        tmp3=q*p';
        if(tmp3>0)
            t=(-0.534*R-e*p')/tmp3;
            tmp4=e+t*q+0.534*R*p;
            counter=counter+(norm(tmp4)<=0.5);
        end  
    end
end
rate4=counter/length(r);


function S=area1(A,B,C,u)                          %һ������¾���ķ����������Ч��� 
%u ������ABC���
%S(1)�������  S(2)��Ч���
global R p
f=0.466*R;     %����
s=-0.3361;
F=f-s;
H=150^2/4/F;                                     %300�ھ�������������߶�
n=50;
o=yuanxin(A,B,C);
w=(B-A)/n; z=(C-A)/n; 
a=[];
for i=1:n
    tmp1=A+(i-1)*w+(w+z)/3;
    tmp2=A+(i-1)*w+2*(w+z)/3;
    a=[a;tmp1+(0:n-i)'*z;tmp2+(0:n-i-1)'*z];            %�ʷ���������
end
tmp=sqrt(sum((a-o).*(a-o),2));
e=o+R./tmp.*(a-o);                                    %�ʷ����������潻�������
s2=e*p'<=H-R+s;                                       %�ж�e�Ƿ��ڹ�������
q=2*(o-e)*p'.*(o-e)/R^2-p;
tmp=q*p';
k=find(tmp>0);
t=zeros(n^2,1);
t(k)=(-0.534*R-e*p')./tmp;
tmp=e(k,:)+t(k).*q(k,:)+0.534*R*p;
s1=sum(tmp.*tmp,2)<=0.25;
S=[sum(s1),sum(s2)]*u/n/n;
end









function S=area2(A,B,C,u)                          %��׼����ķ����������Ч��� 
%u ������ABC���
%S(1)�������  S(2)��Ч���
global R
p=[0,0,1];
H=R-sqrt(R^2-150^2);
n=50;
o=yuanxin(A,B,C);
w=(B-A)/n; z=(C-A)/n;
a=[];
 for i=1:n
    tmp1=A+(i-1)*w+(w+z)/3;
    tmp2=A+(i-1)*w+2*(w+z)/3;
    a=[a;tmp1+(0:n-i)'*z;tmp2+(0:n-i-1)'*z];            %�ʷ���������
end
tmp=sqrt(sum((a-o).*(a-o),2));
e=o+R./tmp.*(a-o);                                    %�ʷ����������潻�������
s2=e(:,3)<=H-R;                                       %�ж�e�Ƿ��ڹ�������
q=2*(o-e)*p'.*(o-e)./R^2-p;
tmp=q*p';
k=find(tmp>0);
t=zeros(n^2,1);
t(k)=(-0.534*R-e*p')./tmp;
tmp=e(k,:)+t(k).*q(k,:)+0.534*R*p;
s1=sum(tmp.*tmp,2)<=0.25;
S=[sum(s1),sum(s2)]*u/n/n;
end







function o=yuanxin(A,B,C)                          %��A��B��C�� Ϊ����ľ���Բ��
R=300.4;
o=zeros(1,3);
p1=A-C;
p2=A-B;
d1=(dot(A,A)-dot(C,C))/2;
d2=(dot(A,A)-dot(B,B))/2;
n=cross(p1,p2);
dc=(d1*p2(3)-d2*p1(3));
da=(d1*p2(1)-d2*p1(1));
db=(d1*p2(2)-d2*p1(2));
if n(1)~=0
    o(1)=0;
    o(2)=dc/n(1);
    o(3)=-db/n(1);
elseif n(2)~=0
    o(2)=0;
    o(1)=-dc/n(2);
    o(3)=da/n(2);
else
    o(3)=0;
    o(1)=db/n(3);
    o(2)=-da/n(3);
end
    a=dot(n,n);
    b=2*dot(n,o-A);
    c=dot(o-A,o-A)-R^2;
    if n(3)>0
        t=(-b+sqrt(b^2-4*a*c))/2/a;
    else
        t=(-b-sqrt(b^2-4*a*c))/2/a;
    end
    o=n*t+o;
end




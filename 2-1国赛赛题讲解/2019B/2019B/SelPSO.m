function [xm, fv]=SelPSO(fitness, N , w,c1, c2,xmax,xmin,M,D)
%������Ȼѡ�������Ⱥ�㷨
% fitness:���Ż���Ŀ�꺯��
% N:������Ŀ
% w: ����Ȩ��
% c1: ѧϰ����1
% c2: ѧϰ����2
% �Ա�������������ֵ��xmax, ������
% �Ա������������Сֵ��xmin  ������
% M: ����������
% D:�Ա�������
% xm�� Ŀ�꺯��ȡ��Сֵʱ���Ա���ֵ
% fv:  Ŀ�꺯����Сֵ
format long
Vm=0.2*(xmax-xmin);
x=zeros(N,D);v=x;
for i=1:N
    for j=1:D
        x(i,j)=xmin(j)+rand*(xmax(j)-xmin(j));
        v(i,j)=Vm(j)*(-1+2*rand);
    end
end
    p=zeros(N,1);y=zeros(N,D);
for i=1:N
    p(i)=fitness(x(i,:));
    y(i,:)=x(i,:);
end
px=x(N,:); pg=p(N);
for i=1:(N-1)
    if p(i)<pg
       pg=p(i);
       px=x(i,:);
    end
end
f=zeros(N,1);
for t=1:M
    for i=1:N
        v(i,:)=w*v(i,:)+c1*(y(i,:)-x(i,:))*rand+c2*rand*(px-x(i,:));
        for j=1:D
            if v(i,j)>Vm(j)
                v(i,j)=Vm(j);
            end
            if v(i,j)<-Vm(j)
                v(i,j)=-Vm(j);
            end
        end
        x(i,:)=x(i,:)+v(i,:);
        tp1=x(i,:)<=xmax;
        tp2=x(i,:)>=xmin;
        if (sum(tp1)+sum(tp2))==2*D
            f(i)=fitness(x(i,:));
        else
            f(i)=inf;
        end
        if f(i)<p(i)
           p(i)=f(i);
           y(i,:)=x(i,:);
        end
        if p(i)<pg
           pg=p(i);
           px=y(i,:);
        end
    end
    [~,sortx]=sort(f);
    exindex=round((N-1)/2);
    x(sortx((N-exindex+1):N))=x(sortx(1:exindex));
   v(sortx((N-exindex+1):N))=v(sortx(1:exindex));
end
xm=px';
fv=pg;
                
            
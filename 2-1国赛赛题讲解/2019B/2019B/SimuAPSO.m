function [xm, fv]=SimuAPSO(fitness, N , w,c1, c2,xmax,xmin,M,D)
%����ģ���˻������Ⱥ�㷨
% fitness:���Ż���Ŀ�꺯��
% N:������Ŀ
% w: �˻���
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
p=zeros(N,1); y=x;
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

tfit=zeros(N,1); f=tfit;comfit=tfit;
T=pg/log(5);
for t=1:M
    groupfit=pg;
    for i=1:N
       tfit(i)=exp(-(p(i)-groupfit)/T);
    end
    sumtfit=sum(tfit);
    tfit=tfit/sumtfit;
    pbet=rand();
    for i=1:N
        comfit(i)=sum(tfit(1:i));
        if pbet<=comfit(i)
            pg_plus=x(i,:);
            break;
        end
    end
    C=c1+c2;
    ksi=2/abs(2-C-sqrt(C^2-4*C));
    for i=1:N
        v(i,:)=ksi*v(i,:)+c1*(y(i,:)-x(i,:))*rand+c2*rand*(pg_plus-x(i,:));
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
    T=T*w;
end
xm=px';
fv=pg;
                
            
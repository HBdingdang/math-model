function [xm, fv]=BreedPSO(fitness, N ,w, c1, c2,xmax,xmin,Pc,Sp,M,D)
%�����ӽ�������Ⱥ�㷨
% fitness:���Ż���Ŀ�꺯��
% N:������Ŀ
% w: ����Ȩ��
% c1: ѧϰ����1
% c2: ѧϰ����2
% �Ա�������������ֵ��xmax, ������
% �Ա������������Сֵ��xmin  ������
% Pc: �ӽ�����
% SP: �ӽ��صĴ�С����
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
f=p;childx1=x;childv1=v;
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
        r1=rand;
        if r1<Pc
           numPool=round(Sp*N);
           PoolX=x(1:numPool,:);
           PoolVx=v(1:numPool,:);
           for j=1:numPool
               seed1=floor(rand*(numPool-1))+1;
               seed2=floor(rand*(numPool-1))+1;
               pb=rand;
               childx1(j,:)=pb*PoolX(seed1,:)+(1-pb)*PoolX(seed2,:);
               childv1(j,:)=(PoolVx(seed1,:)+PoolVx(seed2,:))*norm(PoolVx(seed1,:))/norm(PoolVx(seed1,:)+PoolVx(seed2,:));
           end
           x(1:numPool,:)=childx1(1:numPool,:);
           v(1:numPool,:)=childv1(1:numPool,:);
        end
    end        
end
xm=px';
fv=pg;
                
            
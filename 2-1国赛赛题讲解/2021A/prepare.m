global P e L l R  index side_ind_square ind_area
R=300.4;       %����뾶
f=0.466*R;     %����
s=-0.3361;
F=f-s;
alpha=36.795*pi/180;
beta=78.169*pi/180;
p=[cos(alpha)*cos(beta),cos(beta)*sin(alpha),sin(beta)];
H1=150^2/4/F*1.01^2;                    %�ھ���1%��������߶�
H0=4*(sqrt(F)*(H1+F)^(3/2)-F^2)/3/R;    %�����α������߶�
filename="����1.csv";
[E,c]=xlsread(filename,1);             %E:�����ڵ�ԭʼ����
ind2=string(c(2:length(c),1));          %ind2  �����ڵ��ʶ
k=find(E(:,1:3)*p'<=H0-R);              %�ڶ����α������ڵĽڵ�����
filename="����2.csv";
G=xlsread(filename,1);
X=E-G(:,4:6);            
l=sum(X.*X,2);                          %����������ƽ��
X=G(:,4:6)-G(:,1:3);
L=sqrt(sum(X.*X,2));                    %�ٶ���ԭ��      
e=(G(:,4:6)-G(:,1:3))./L;               %�ٶ�������
P=G(:,1:3);                             %�ٶ����¶�����
filename="����3.csv";
[~,c]=xlsread(filename,1);
tmp=string(c(2:length(c),1:3));
ind_area=zeros(4300,4);                   %���涥���������������
ind=zeros(4300,3);                        %���涥����������
for i=1:2226                            
    j=find(tmp==ind2(i));
    ind(j)=i;
end
ind_area(:,1:3)=ind;
tmp=zeros(2226);
for i=1:4300
    tmp(ind(i,1),ind(i,2))=1;
    tmp(ind(i,2),ind(i,3))=1;
    tmp(ind(i,3),ind(i,1))=1;
end
for i=1:2226
    for j=i+1:2226
        if tmp(j,i)
            tmp(i,j)=1;
            tmp(j,i)=0;
        end 
    end
end
t=find(tmp);  
n=length(t);                        %n=6525  ��������
side_ind_square=zeros(n,3);         %1,2 �������˵㡢3 ��������
for i=1:n
    a=floor(t(i)/2226)+1;
    b=mod(t(i),2226);
    side_ind_square(i,1:2)=[b,a];
    side_ind_square(i,3)=sum((E(b,:)-E(a,:)).^2);
end

index=zeros(2226,n);                %|Q_i-Q_j|^2���ݶ�ʱ����ľ���             
for i=1:n
    index(side_ind_square(i,1),i)=1;
    index(side_ind_square(i,2),i)=-1;
end


for i=1:4300
     tmp=find(side_ind_square(:,1)==ind(i,1));
     tmp1= side_ind_square(tmp,2)==ind(i,2);
     a=sqrt(side_ind_square(tmp(tmp1),3));
     tmp1= side_ind_square(tmp,2)==ind(i,3);
     b=sqrt(side_ind_square(tmp(tmp1),3));
     tmp=find(side_ind_square(:,1)==ind(i,2));
     tmp1=find(side_ind_square(tmp,2)==ind(i,3));
     c=sqrt(side_ind_square(tmp(tmp1),3));
     p=(a+b+c)/2;
     ind_area(i,4)=sqrt(p*(p-a)*(p-b)*(p-c));
end

save data.mat E L P e l ind_area ind2 index side_ind_square;



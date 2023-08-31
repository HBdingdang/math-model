global R E
load('data.mat','l','L','P','E','e');

R=300.4;       %����뾶
%l                    ����������ƽ��
%L                    �ٶ���ԭ��  
%e                    �ٶ�������
%P                    �ٶ����¶�����
%E                    �����ڵ�ԭ����

figure
ax1=subplot(1,3,1);
ax2=subplot(1,3,2);
ax3=subplot(1,3,3);
s=-0.6:0.1:0.6; S=s;
for i=1:length(s)
    S(i)=shensuoliang(s(i));
end
plot(ax1,s,S);
ylabel(ax1,'�α������дٶ�����������������ף�');
[s1,v1]=fminbnd(@shensuoliang, -0.6,0.6);   


s=-0.6:0.1:0.6; S=s;
for i=1:length(s)
    S(i)=weiyiliang(s(i));
end
plot(ax2,s,S);
ylabel(ax2, '�α������������ƶ��������ƽ������^2��');
xlabel(ax2,'x�᣺��͵�Ĵٶ��������������ף�');
[s2,v2]=fminbnd(@weiyiliang, -0.6,0.6); 

for i=1:length(s)
    S(i)=pingjunjuli(s(i));
end
plot(ax3,s,S);
ylabel(ax3,'�α���������Ӧ���ƶ���ƽ�����루�ף�');


[s3,v3]=fminbnd(@pingjunjuli, -0.6,0.6);   


function S=pingjunjuli(s)          %s ��͵�Ĵٶ����쳤��
global R 
f=0.466*R;     %����
F=f-s;
h=0.01;
xp=0:h:150;
hp=xp.^2/4/F;
zp=s-R+hp;
hb=4*(sqrt(F)*(hp+F).^1.5-F^2)/3/R;
zb=hb-R;
xb=sqrt(R^2-zb.^2);
tmp=sqrt((xp-xb).^2+(zp-zb).^2);
S=2*sum(tmp)*h/150;
end

function S=shensuoliang(s)           %s ��͵�Ĵٶ����쳤��
global R E L l e P 
f=0.466*R;     %����
F=f-s;
H1=150^2/4/F;                       %�ھ�300��������߶�
H0=4*(sqrt(F)*(H1+F)^(3/2)-F^2)/3/R; %�����α������߶�
k=find(E(:,3)<=H0-R);                %�α������ڵĽڵ�����
X=E(k,:);                            %�α������ڵĽڵ�ԭ����
t=((X(:,3)+R)~=0).*height(F,X(:,3)+R);
tmp=sqrt(4*F*t);
t1=X(:,1:2);
t2=sqrt(sum(t1.*t1,2));
tmp1=find(t2);
t1(tmp1,:)=t1(tmp1,:)./t2(tmp1);
Q=[tmp.*t1,t+s-R];                  %�α������ڵĽڵ�������
tmp2=P(k,:)-Q;
b=sum(e(k,:).*tmp2,2);
c=sum(tmp2.*tmp2,2)-l(k);
tmp3=sqrt(b.^2-c);
S=max(min(abs(-L(k)-b-tmp3), abs(-L(k)-b+tmp3)));  %�α�����ٶ��������������ֵ
end 

function S=weiyiliang(s)           %s ��͵�Ĵٶ����쳤��
global R E
f=0.466*R;     %����
F=f-s;
H1=150^2/4/F;                        %�ھ�300��������߶�
H0=4*(sqrt(F)*(H1+F)^(3/2)-F^2)/3/R; %�����α������߶�
k=find(E(:,3)<=H0-R);                %�α������ڵĽڵ�����
X=E(k,:);                            %�α������ڵĽڵ�ԭ����
t=((X(:,3)+R)~=0).*height(F,X(:,3)+R);
tmp=sqrt(4*F*t);
t1=X(:,1:2);
t2=sqrt(sum(t1.*t1,2));
tmp1=find(t2);
t1(tmp1,:)=t1(tmp1,:)./t2(tmp1);
Q=[tmp.*t1,t+s-R];                  %�α������ڵĽڵ�������
tmp=Q-X;
S=max(sum(tmp.*tmp,2));
end 


function z=height(f,x)             %��x�����������α�Ϊ��Ϊz,����Ϊf��������
global R 
z=((3.*x*R/4+f^2)./sqrt(f)).^(2/3)-f;
end




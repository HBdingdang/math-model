global R E
load('data.mat','l','L','P','E','e');

R=300.4;       %球面半径
%l                    下拉索长度平方
%L                    促动器原长  
%e                    促动器方向
%P                    促动器下端坐标
%E                    主索节点原坐标

figure
ax1=subplot(1,3,1);
ax2=subplot(1,3,2);
ax3=subplot(1,3,3);
s=-0.6:0.1:0.6; S=s;
for i=1:length(s)
    S(i)=shensuoliang(s(i));
end
plot(ax1,s,S);
ylabel(ax1,'形变区域中促动器的最大伸缩量（米）');
[s1,v1]=fminbnd(@shensuoliang, -0.6,0.6);   


s=-0.6:0.1:0.6; S=s;
for i=1:length(s)
    S(i)=weiyiliang(s(i));
end
plot(ax2,s,S);
ylabel(ax2, '形变区域中索点移动的最大量平方（米^2）');
xlabel(ax2,'x轴：最低点的促动器的伸缩量（米）');
[s2,v2]=fminbnd(@weiyiliang, -0.6,0.6); 

for i=1:length(s)
    S(i)=pingjunjuli(s(i));
end
plot(ax3,s,S);
ylabel(ax3,'形变区域中相应点移动的平均距离（米）');


[s3,v3]=fminbnd(@pingjunjuli, -0.6,0.6);   


function S=pingjunjuli(s)          %s 最低点的促动器伸长量
global R 
f=0.466*R;     %焦距
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

function S=shensuoliang(s)           %s 最低点的促动器伸长量
global R E L l e P 
f=0.466*R;     %焦距
F=f-s;
H1=150^2/4/F;                       %口径300的抛物面高度
H0=4*(sqrt(F)*(H1+F)^(3/2)-F^2)/3/R; %参与形变的球面高度
k=find(E(:,3)<=H0-R);                %形变区域内的节点索引
X=E(k,:);                            %形变区域内的节点原坐标
t=((X(:,3)+R)~=0).*height(F,X(:,3)+R);
tmp=sqrt(4*F*t);
t1=X(:,1:2);
t2=sqrt(sum(t1.*t1,2));
tmp1=find(t2);
t1(tmp1,:)=t1(tmp1,:)./t2(tmp1);
Q=[tmp.*t1,t+s-R];                  %形变区域内的节点新坐标
tmp2=P(k,:)-Q;
b=sum(e(k,:).*tmp2,2);
c=sum(tmp2.*tmp2,2)-l(k);
tmp3=sqrt(b.^2-c);
S=max(min(abs(-L(k)-b-tmp3), abs(-L(k)-b+tmp3)));  %形变区域促动器伸缩量的最大值
end 

function S=weiyiliang(s)           %s 最低点的促动器伸长量
global R E
f=0.466*R;     %焦距
F=f-s;
H1=150^2/4/F;                        %口径300的抛物面高度
H0=4*(sqrt(F)*(H1+F)^(3/2)-F^2)/3/R; %参与形变的球面高度
k=find(E(:,3)<=H0-R);                %形变区域内的节点索引
X=E(k,:);                            %形变区域内的节点原坐标
t=((X(:,3)+R)~=0).*height(F,X(:,3)+R);
tmp=sqrt(4*F*t);
t1=X(:,1:2);
t2=sqrt(sum(t1.*t1,2));
tmp1=find(t2);
t1(tmp1,:)=t1(tmp1,:)./t2(tmp1);
Q=[tmp.*t1,t+s-R];                  %形变区域内的节点新坐标
tmp=Q-X;
S=max(sum(tmp.*tmp,2));
end 


function z=height(f,x)             %高x的球面等面积形变为高为z,焦距为f的抛物面
global R 
z=((3.*x*R/4+f^2)./sqrt(f)).^(2/3)-f;
end




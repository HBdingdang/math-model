function [t,theta1,theta]=drum_angle(F,h,T)
% F: ����������ʱ���б�
% h: �����½��߶�
% T: ����ʱ��
% t: ʱ��ڵ�����
% theta1: �����������
% theta: t=0.1ʱ��������б��
m=8;L=1.7; r0=0.2; dt=0.001; C=1;rho=1.2258;g=9.8;M=3.6;H=0.11;
F0=M*g*L/(m*h);    % ����ƽ��λ��ʱ���ӵ�����
S=pi*r0^2;   % �������
I=M*(r0^2/2+H^2/3);  % ����x���y���ת������I
J=M*r0^2;               % ����z���ת������J
X=zeros(m,3);   %X ������꣬ ��һ�������x������
x=zeros(m,3); % x ��Ӧ�����Ķ�����
r1=zeros(m,3);  % ���ĵ����Ķ˵�����

t=min([0,F(1,:)]):dt:T;
e1=zeros(length(t),3); e2=e1;e3=e1;     %����x,y,z�ᵥλ�����ڴ�
w=e1;                                   %������ٶ������ڴ�
e1(1,:)=[1,0,0]; e2(1,:)=[0,1,0]; e3(1,:)=[0,0,1];   %x,y,z���ʼ��λ����
v=zeros(length(t),3);                   % �����ƽ���ٶ������ڴ�        
c=v;                                    %������������ڴ�;
f=zeros(m,3);
for i=1:length(t)-1
    
    for j=1:m
        r1(j,:)=r0*(e1(i,:)*cos(2*pi*(j-1)/m)+e2(i,:)*sin(2*pi*(j-1)/m));
    end
    
    
    for j=1:m
        x(j,:)=r1(j,:)+c(i,:);
    end
    
    for j=1:m
      tmp=dot([x(j,1),x(j,2)],[cos(2*pi*(j-1)/m),sin(2*pi*(j-1)/m)])+sqrt(L^2-x(j,1)^2-x(j,2)^2-(x(j,3)-h)^2);  
      X(j,:)=[tmp*[cos(2*pi*(j-1)/m),sin(2*pi*(j-1)/m)],h];
    end
    
    
    Q=F0+(F(1,:)<=t(i)).*(F(2,:)-F0);
    for j=1:m
        f(j,:)=(Q(j))*(X(j,:)-x(j,:))/norm(X(j,:)-x(j,:));    %������
    end
    
    
    tmp1=sum(f)/M-[0,0,g]-C*rho*S*norm(v(i,:))*v(i,:)/(2*M);
    vp=v(i,:)+tmp1*dt;                                   %�м��ٶ����������
    cp=c(i,:)+v(i,:)*dt;
     
       
    M0=zeros(1,3);   %��������
    for j=1:8
        M0=M0+cross(r1(j,:),f(j,:));
    end
    
    M1=dot(M0,e1(i,:));                   % M��e1 ����ķ���
    M2=dot(M0,e2(i,:));                   % M��e2 ����ķ���
    M3=dot(M0,e3(i,:));                   % M��e3 ����ķ���
    
    
    tmp=-(J/I-1)*w(i,2)*w(i,3)+M1/I;
    w1p=w(i,1)+tmp*dt;
    
    tmp=(J/I-1)*w(i,1)*w(i,3)+M2/I;
    w2p=w(i,2)+tmp*dt;
    
    tmp=M3/J;
    w3p=w(i,3)+tmp*dt;                           %�½��ٶ������� w(1)*e1+w(2)*e2+w3*e3
    
    tmp=cross(w(i,2)*e2(i,:)+w(i,3)*e3(i,:),e1(i,:));
    e1p=(e1(i,:)+tmp*dt)/(norm(e1(i,:)+tmp*dt));
    
    tmp=cross(w(i,1)*e1(i,:)+w(i,3)*e3(i,:),e2(i,:));
    e2p=(e2(i,:)+tmp*dt)/(norm(e2(i,:)+tmp*dt));
    
    tmp=cross(w(i,1)*e1(i,:)+w(i,2)*e2(i,:),e3(i,:));
    e3p=(e3(i,:)+tmp*dt)/(norm(e3(i,:)+tmp*dt));     %�µ�x,y,z��
    
    for j=1:m
        r1(j,:)=r0*(e1p*cos(2*pi*(j-1)/m)+e2p*sin(2*pi*(j-1)/m));
    end
    
    for j=1:m
        x(j,:)=r1(j,:)+cp;
    end
    
    
    for j=1:m
      tmp=dot([x(j,1),x(j,2)],[cos(2*pi*(j-1)/m),sin(2*pi*(j-1)/m)])+sqrt(L^2-x(j,1)^2-x(j,2)^2-(x(j,3)-h)^2);  
      X(j,:)=[tmp*[cos(2*pi*(j-1)/m),sin(2*pi*(j-1)/m)],h];
    end                                                             %�µ��ֶ�����
    
    for j=1:m
        f(j,:)=(Q(j))*(X(j,:)-x(j,:))/norm(X(j,:)-x(j,:));    
    end
    tmp1=sum(f)/M-[0,0,g]-C*rho*S*norm(vp)*vp/(2*M);
    vc=v(i,:)+tmp1*dt;                                   %�м��ٶ����������
    cc=c(i,:)+vp*dt;
    
     M0=zeros(1,3);   %��������
    for j=1:m
        M0=M0+cross(r1(j,:),f(j,:));
    end
    
    M1=dot(M0,e1p);                   % M��e1 ����ķ���
    M2=dot(M0,e2p);                   % M��e2 ����ķ���
    M3=dot(M0,e3p);                   % M��e3 ����ķ���
    
    
    tmp=-(J/I-1)*w2p*w3p+M1/I;
    w1c=w(i,1)+tmp*dt;
    
    tmp=(J/I-1)*w1p*w3p+M2/I;
    w2c=w(i,2)+tmp*dt;
    
    tmp=M3/J;
    w3c=w(i,3)+tmp*dt;                           %�½��ٶ������� w(1)*e1+w(2)*e2+w3*e3
    
    tmp=cross(w2p*e2p+w3p*e3p,e1p);
    e1c=(e1(i,:)+tmp*dt)/(norm(e1(i,:)+tmp*dt));
    
    tmp=cross(w1p*e1p+w3p*e3p,e2p);
    e2c=(e2(i,:)+tmp*dt)/(norm(e2(i,:)+tmp*dt));
    
    tmp=cross(w1p*e1p+w2p*e2p,e3p);
    e3c=(e3(i,:)+tmp*dt)/(norm(e3(i,:)+tmp*dt));     %�µ�x,y,z��
    
    v(i+1,:)=(vp+vc)/2;
    c(i+1,:)=(cc+cp)/2;
    w(i+1,1)=(w1p+w1c)/2;
    w(i+1,2)=(w2p+w2c)/2;
    w(i+1,3)=(w3p+w3c)/2;
    e1(i+1,:)=(e1p+e1c)/2;
    e2(i+1,:)=(e2p+e2c)/2;
    e3(i+1,:)=(e3p+e3c)/2;
    
    
end
theta1=zeros(length(t),1);    
for i=1:length(t)
    theta1(i)=acos(e3(i,3))*180/pi;
end
theta=theta1(find(t==0.1));



end
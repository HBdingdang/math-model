function speed = min_descend_speed( a,b,v0,e)
%ascend�� ����������
%  a: �������С�߶�
%  b: ����ײ�߶�
%  v0: ������ײλ�õ��ٶ�
% speed: ʹ���������b����ײ������ĸ߶ȳ���0.4�׵���С�½��ٶ�
% M: ͬ�Ĺĵ�����
% m: ������
% e: ��ײ�ָ�ϵ��  ��ȡe=1
M=3.6; m=0.27;
[~,~,v1]=descend(a,b);
v=min_ascend_speed(a,b);
tmp=((m+M)*v-M*(1+e)*v0)/(m-e*M);
speed=max(-1*(tmp<0)*tmp,abs(v1));
end


clc
clear all
close all

%%
N=200000;%�źų���
M=4;%���ƽ���
d = randi([0 M-1],N,1);
y = qammod(d,M,'PlotConstellation',true);
y=y.';%����IQ
snr=20;
yn=awgn(y,snr,'measured');

%% channel
h=[1,0.2,0.3];%�ŵ�����,�ŵ�Խ���ӣ���Ҫ�����˲�������
yh=conv(yn,h)/sum(h);
yh=yh(length(h):end);%�����ŵ�����ź�
%% CMA
%parameter
yh=yh/max(abs(yh))*2;
u=yh;
R2=mean(abs(u).^4)/mean(abs(u).^2);


L=21;%�˲�������
mu=0.001;%����
U=zeros(L,N);
u=[zeros(1,L-1),u];

for ii=1:N
    U(:,ii)=u(ii:ii+L-1);
end


%init0
w=zeros(1,L)';%�˲���ϵ��
w(round(L/2))=1;
e=zeros(1,N);
s=zeros(1,N);
d=zeros(1,N);
%get_g()��CMA�ķ����Ժ���
e(1)=get_g(s(1),R2,2)-w'*U(:,1);
for ii=2:N
    s(ii)=w'*U(:,ii);
    d(ii)=get_g(s(ii),R2,2);
    e(ii)=d(ii)-s(ii);
    w=w+mu*U(:,ii)*conj(e(ii));
end


figure
subplot(2,2,1)
plot(yn,'x')
title('ֻ�и�˹����������')
subplot(2,2,2)
plot(yh,'x')
title('�����ŵ��������')
subplot(2,2,3)
plot(s(end-1000:end),'x')
title('CMA�����������')
subplot(2,2,4)
plot(abs(e))
title('��������')

%�źŹ�һ����Ա�
s=s/mean((abs(s)));
yn=yn/mean((abs(yn)));

%ȷ���ź��ӳ�
e_mean=zeros(1,L);
for ii=1:L
    e_mean(ii)=mean(abs(s(1+ii-1:end)-yn(1:end-ii+1)));
end
[~,min_pos]=min(e_mean);

figure
plot(abs(s(1+min_pos-1:end)-yn(1:end-min_pos+1)))
title('�������')
disp(["����ź���Ҫ�ӳ�" ,min_pos, "����λ������ԭ�ź���ͬ"])

figure
subplot(211)
plot(yn,'x');
subplot(212)
plot(s(end-1000:end),'x')
suptitle('ԭʼ�źźͻָ��źŶԱ�')
% e0=abs(s-yn);
% e1=abs(s(2:end)-yn(1:end-1));
% e2=abs(s(3:end)-yn(1:end-2));
% e3=abs(s(4:end)-yn(1:end-3));
% 
% % �����һ��������߻�������˵���㷨��Ч
% figure
% subplot(221)
% plot(e0)
% title('����λ���')
% subplot(222)
% plot(e1)
% title('��λ1���������')
% subplot(223)
% plot(e2)
% title('��λ2���������')
% subplot(224)
% plot(e3)
% title('��λ3���������')
% suptitle('4����λ���')

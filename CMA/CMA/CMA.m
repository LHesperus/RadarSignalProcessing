clc
clear all
close all

%%
N=200000;%信号长度
M=4;%调制阶数
d = randi([0 M-1],N,1);
y = qammod(d,M,'PlotConstellation',true);
y=y.';%基带IQ
snr=20;
yn=awgn(y,snr,'measured');

%% channel
h=[1,0.2,0.3];%信道参数,信道越复杂，需要增加滤波器阶数
yh=conv(yn,h)/sum(h);
yh=yh(length(h):end);%经过信道后的信号
%% CMA
%parameter
yh=yh/max(abs(yh))*2;
u=yh;
R2=mean(abs(u).^4)/mean(abs(u).^2);


L=21;%滤波器长度
mu=0.001;%步长
U=zeros(L,N);
u=[zeros(1,L-1),u];

for ii=1:N
    U(:,ii)=u(ii:ii+L-1);
end


%init0
w=zeros(1,L)';%滤波器系数
w(round(L/2))=1;
e=zeros(1,N);
s=zeros(1,N);
d=zeros(1,N);
%get_g()是CMA的非线性函数
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
title('只有高斯噪声的星座')
subplot(2,2,2)
plot(yh,'x')
title('经过信道后的星座')
subplot(2,2,3)
plot(s(end-1000:end),'x')
title('CMA收敛后的星座')
subplot(2,2,4)
plot(abs(e))
title('收敛曲线')

%信号归一化后对比
s=s/mean((abs(s)));
yn=yn/mean((abs(yn)));

%确定信号延迟
e_mean=zeros(1,L);
for ii=1:L
    e_mean(ii)=mean(abs(s(1+ii-1:end)-yn(1:end-ii+1)));
end
[~,min_pos]=min(e_mean);

figure
plot(abs(s(1+min_pos-1:end)-yn(1:end-min_pos+1)))
title('误差曲线')
disp(["输出信号需要延迟" ,min_pos, "个单位才能与原信号相同"])

figure
subplot(211)
plot(yn,'x');
subplot(212)
plot(s(end-1000:end),'x')
suptitle('原始信号和恢复信号对比')
% e0=abs(s-yn);
% e1=abs(s(2:end)-yn(1:end-1));
% e2=abs(s(3:end)-yn(1:end-2));
% e3=abs(s(4:end)-yn(1:end-3));
% 
% % 如果有一个误差曲线会收敛，说明算法有效
% figure
% subplot(221)
% plot(e0)
% title('不移位误差')
% subplot(222)
% plot(e1)
% title('移位1个数据误差')
% subplot(223)
% plot(e2)
% title('移位2个数据误差')
% subplot(224)
% plot(e3)
% title('移位3个数据误差')
% suptitle('4种移位误差')

%%%%%%%%%%%%%% ADMM for Lasso %%%%%%%%%%%%%%%%%
clc
clear
close all
% simulated signal
SampFreq = 2000;
t = 0:1/SampFreq:1;
a1 = sawtooth(2*pi*20*t,0.1)+2;
Sig = a1.*cos(2*pi*(0.8+50*t + 525*t.^2 -300*t.^3));
IF = 50 + 1050*t - 900*t.^2;

%% Original signal
fsize = 12;
figure;
subplot(3,1,1)
plot(t,Sig);
xlabel({'Time (s)','(a)'});
ylabel('Amplitude');
ylim([-4 4]);
set(gca,'fontsize',fsize,'linewidth',1)

subplot(3,1,2)
plot(t,a1);
xlabel({'Time (s)','(b)'});
ylabel('Amplitude');
ylim([0 4]);
legend('IA');
set(gca,'fontsize',fsize,'linewidth',1)

subplot(3,1,3)
plot(t,IF);
xlabel({'Time (s)','(c)'});
ylabel('Frequency (Hz)');
ylim([0 400]);
legend('IF');
set(gca,'fontsize',fsize,'linewidth',1)

%% Time-frequency representation

%输入参数
stft_type ='M';
tf_type ='SET1';
direction = 'F';
gamma = 0.001;
s = 0.02;
tradeoff = 0.009;
delta = 60;

%计算TFR
[Wx,TFx,Ifreq,GD,Rep,Chirprate,q,t,f] = TFM(Sig,SampFreq,s,stft_type,tf_type,gamma);
%TFR绘制
figure
% subplot(2,1,1)
imagesc(t,f,abs(Wx'));
axis xy
xlabel('Time (s)','FontSize',12);
ylabel('Frequency (Hz)','FontSize',12);
ylim([0 400]);
title('STFT');
set(gca,'FontSize',12);
% subplot(2,1,2)
% imagesc(t,f,abs(TFx'));
% axis xy
% xlabel('Time (s)','FontSize',20);
% ylabel('Frequency (Hz)','FontSize',20);
% title(['stft-type:',stft_type,'      tf-type:',tf_type]);
% set(gca,'FontSize',20);
%% parameter setting
iniIF = 300*ones(1,length(t));% initial guess for the IFs for the three signal modes
beta = 1e-6; % this parameter can be smaller which will be helpful for the convergence, but it may cannot properly track fast varying IFs
tol = 1e-8;%
lambda =5e-4;
% tic
[IFmset1,IA1,smset1] = NCME(Sig,SampFreq,iniIF,lambda,beta,tol);
% toc
 
    %% estimated IF,IA,S
    figure;
    subplot(3,1,1)
    plot(t,smset1(1,:,end),'r','linewidth',1)	% estimated mode
    hold on;
    plot(t,Sig,'b','linewidth',1);
    xlabel({'Time/s','(a)'});
    ylabel('Amplitude');
    ylim([-4 4]);
   set(gca,'fontsize',fsize,'linewidth',1)
   legend('Estimated','Original');
    
    subplot(3,1,2)
    plot(t,IFmset1(:,:,end),'r','linewidth',1); % finally estimated IFs
    hold on
    plot(t,IF,'b','linewidth',1);
    xlabel({'Time/s','(b)'});
    ylabel('Frequency/Hz');
    ylim([0 400]);
    set(gca,'fontsize',fsize,'linewidth',1);
    legend('Estimated','Original');
       
       x1=0.2;
       y1=0.8;
       x2=0.26;
       y2=3.2;
     subplot(3,1,3)
     plot(t,IA1,'r','linewidth',1);
     hold on;
     plot(t,a1,'b','linewidth',1); 
    xlabel({'Time/s','(c)'});
    ylabel('Amplitude');
    ylim([0 4]);
    rectangle('Position',[x1 y1 x2-x1 y2-y1],'EdgeColor','r','Linewidth',1);
    legend('Estimated','Original');
    set(gca,'fontsize',fsize,'linewidth',1)
    
%     h1=axes('position',[0.7 0.15 0.2 0.1]);%让 坐标轴的左下角 与 窗口左侧 的距离时窗口宽度的8%，距离下侧10%
%                                            %整个坐标轴的宽占85%，高占85%。一个小框就出来了
%      axis(h1);
%      plot(t,IA1,'r','linewidth',1);
%      hold on;
%      plot(t,a1,'b','linewidth',1); 
%      xlim([x1 x2]);ylim([y1 y2]);
%      set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
    %%  estimated amplitude 估计幅值
     figure;
%    set(gcf,'Position',[50 50 560 840]);
    plot(t,a1,'linewidth',1); 
     hold on;
     plot(t,IA1,'linewidth',1);
    xlabel({'Time/s','(a)'});
    ylabel('Amplitude');
    ylim([0 4]);
    rectangle('Position',[x1 y1 x2-x1 y2-y1],'EdgeColor','k','Linewidth',1);
    legend('Original','Estimated');
    set(gca,'fontsize',fsize,'linewidth',1)
    
     h1=axes('position',[0.60 0.2 0.3 0.15]);%让 坐标轴的左下角 与 窗口左侧 的距离时窗口宽度的8%，距离下侧10%
                                             %整个坐标轴的宽占85%，高占85%。一个小框就出来了
     axis(h1);
     plot(t,a1,'linewidth',1.5); 
     hold on;
     plot(t,IA1,'linewidth',1.5);
     xlim([x1 x2]);ylim([y1 y2]);
     set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
     set(gca,'fontsize',fsize,'linewidth',1)
    
    %% Reconstructed signal重构信号
    figure;
    plot(t,Sig,'linewidth',1); 
    hold on;
    plot(t,smset1(1,:,end),'linewidth',1);
    xlabel({'Time/s','(a)'});
    ylabel('Amplitude');
    ylim([-3 3]);
    legend('Original','Estimated');
    set(gca,'fontsize',fsize,'linewidth',1)
    

%% 评估因子
SNR(a1,IA1)
SNR(IF,IFmset1(:,:,end))
SNR(Sig,smset1(1,:,end))

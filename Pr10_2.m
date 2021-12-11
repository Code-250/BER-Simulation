clear all,  close all
SNRdB =   [0:35];
iter=1000; 
C1=ergcap(SNRdB,1,1,iter);
C2=ergcap(SNRdB,2,2,iter);
C3=ergcap(SNRdB,4,4,iter);
figure, plot(SNRdB,C1(1,:),'b-o', SNRdB,C2(1,:),'b-*', SNRdB,C3(1,:),'b-+');
xlabel('SNR(dB)'); ylabel('bps/Hz'); set(gca,'fontsize',10);
legend('{\it N_T}=1,{\it N_R}=1','{\it N_T}=2,{\it N_R}=2', '{\it N_T}=4,{\it N_R}=4')

%Program  to compute the capacity of  iid and correlated MIMO channels
clear all, close all;
SNRdB=[0:5:25];  SNRlin=10.^(SNRdB/10);
iterations=1000;  nT=2;  nR=2; 
n=min(nT,nR);  I = eye(n);  sq2=sqrt(0.5);
R=[1  0.76*exp(0.17j*pi) ; 0.76*exp(-0.17j*pi)   1 ];
Ciid=zeros(1,length(SNRdB));  Ccorr=zeros(1,length(SNRdB));
for iter=1:iterations
   H = sq2*(randn(nR,nT)+j*randn(nR,nT));
   Hc = H*R^(1/2);
   temp1 =  H'*H/nT;  temp2 = Hc'*Hc /nT;
   for i=1:length(SNRdB)
      Ciid(i) = Ciid(i) + log2(det(I+SNRlin(i)*temp1));
    Ccorr(i) =  Ccorr(i) + log2(det(I+SNRlin(i)*temp2));
   end
end
Ciid = real(Ciid)/iterations;  Ccorr = real( Ccorr)/iterations;
plot(SNRdB,Ciid, SNRdB, Ccorr,':');
xlabel('SNR (dB)'); ylabel('Capacity(bps/Hz)'); set(gca,'fontsize',10)
legend('iid 2x2 fading channels','correlated 2x2 fading channels');

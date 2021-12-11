%Program  to compute the capacity of  random MIMO channels with and  without CSI with  uniform  correlated matrices
clear all,  close all
SNRdB =   [0:25]; iter=1000; nT=4; nR=4;
n=min(nT,nR);  I = eye(n);SNRlin=10.^(SNRdB/10.);rho=0.2;
Rtx=[1      rho     rho^2   rho^3;  rho     1      rho    rho^2;
    rho^2   rho     1       rho;   rho^3   rho^2   rho     1];
Rrx=Rtx;
C1(1,:) = zeros(1,length(SNRdB));
C2(1,:) = zeros(1,length(SNRdB));
   for ii=1:iter
  Hw = sqrt(0.5)*(randn(nR,nT)+j*randn(nR,nT)); 
H = Rrx^(1/2)*Hw*Rtx^(1/2);tmp = H'*H/nT;HH= H'*H; Lamda = svd(H'*H);
  for i=1:length(SNRdB) %random channel generation
         C1(1,i) = C1(1,i)+ log2(real(det(I+SNRlin(i)*tmp)));
         gama = WaterFilling(Lamda,SNRlin(i),nT);
 C2(1,i) = C2(1,i)+log2(real(det(I+(SNRlin(i)/nT)*diag( gama )*diag(Lamda))));
      end
   end
C1 = C1/iter; C2=C2/iter;
figure, plot(SNRdB,C1(1,:),'-', SNRdB,C2(1,:),'--');
xlabel('SNR(dB)'); ylabel('bps/Hz'); set(gca,'fontsize',10);
legend('with unknown CSI at the transmitter' ,'with known CSI at the transmit-ter')

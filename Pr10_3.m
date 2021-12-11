%Program to compute the capacity of random MIMO channels without CSI with uniform and non-uniform correlated matrices
clear all,  close all
SNRdB =   [0:25];iter=1000; nT=4; nR=4;
n=min(nT,nR);  I = eye(n);SNRlin=10.^(SNRdB/10.);rho=0.2;
Rtx=[1      0.3169     0.3863    0.0838;  0.3169    1      0.7128   0.5626;
     0.3863   0.7128     1       0.5354;   0.0838   0.5626   0.5354   1];
Rrx=[1      0.1317    0.1992  0.2315;  0.1317     1       0.1493     0.1907;
    0.1992   0.1493     1       0.1996;  0.2315  0.1907  0.1996   1];
Rtxu=[1      rho     rho^2   rho^3;  rho     1      rho    rho^2;
    rho^2   rho     1       rho;   rho^3   rho^2   rho     1];
Rrxu= Rtxu;
C1(1,:) = zeros(1,length(SNRdB)); C2(1,:) = zeros(1,length(SNRdB));
   for ii=1:iter
      Hw = sqrt(0.5)*(randn(nR,nT)+j*randn(nR,nT)); 
      H = Rrx^(1/2)*Hw*Rtx^(1/2);  tmp = H'*H/nT;
      Hu = Rrxu^(1/2)*Hw*Rtxu^(1/2);  tmp1 = Hu'*Hu/nT;;
           for i=1:length(SNRdB) %random channel generation
         C2(1,i) = C2(1,i)+log2(real(det(I+SNRlin(i)*tmp1)));
         C1(1,i) = C1(1,i)+log2(real(det(I+SNRlin(i)*tmp)));
      end
   end
C1 = C1/iter; C2=C2/iter;
figure, plot(SNRdB,C1(1,:),'-', SNRdB,C2(1,:),'--');
xlabel('SNR(dB)'); ylabel('bps/Hz'); set(gca,'fontsize',10);
legend('with non-uniform correlation matrices' ,'with uniform correlation matrices')

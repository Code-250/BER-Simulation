function C=ergcap(SNRdB,nT,nR,iter)
n=min(nT,nR);  I = eye(n);SNRlin=10.^(SNRdB/10.);
   C(1,:) = zeros(1,length(SNRdB));
   for ii=1:iter
      H = sqrt(0.5)*(randn(nR,nT)+j*randn(nR,nT));  
      if nR>=nT,  HH = H'*H;  else  HH = H*H';  end
      for i=1:length(SNRdB) %random channel generation
         C(1,i) = C(1,i)+log2(real(det(I+SNRlin(i)/nT*HH)));
      end
   end
C = C/iter;
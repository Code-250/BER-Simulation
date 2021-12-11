%channel_estimation.m
% for LS/DFT Channel Estimation with linear/spline interpolation
clear all; close all; 
N_fft=32;  N_g=N_fft/8;  N_ofdm=N_fft+N_g;  N_sym=100;
N_ps=4; N_p=N_fft/N_ps; N_d=N_fft-N_p; % Pilot spacing, Numbers of pilots and data per OFDM symbol
N_bps=4; M=2^N_bps; % Number of bits per (modulated) symbol
mod_object = modem.qammod('M',M, 'SymbolOrder','gray');
demod_object = modem.qamdemod('M',M, 'SymbolOrder','gray');
Es=1; A=sqrt(3/2/(M-1)*Es); % Signal energy& QAM normalization factor
EbN0s = [0:5:30];  sq2=sqrt(2);
for i=1:length(EbN0s)
   EbN0 = EbN0s(i); 
   rand('seed',1); randn('seed',1);
   MSE_LSi = 0;  MSE_DFTi=0;
   for nsym=1:N_sym
      X_p = 2*(randn(1,N_p)>0)-1;    % Pilot sequence generation
      msg_int=randint(1,N_fft-N_p,M);    % bit generation
      Data = modulate(mod_object,msg_int)*A;
      ip = 0;    pilot_loc = [];
      for k=1:N_fft
         if mod(k,N_ps)==1
           X(k) = X_p(floor(k/N_ps)+1); pilot_loc = [pilot_loc k]; ip = ip+1;
          else        X(k) = Data(k-ip);
         end
      end
      x = ifft(X,N_fft);     % IFFT
      xt = [x(N_fft-N_g+1:N_fft) x];       % Add CP
      h = [(randn+j*randn) (randn+j*randn)/2]; % generates a (2-tap) channel
      H = fft(h,N_fft); channel_length = length(h); % True channel and its time-domain lengt
      y_channel = conv(xt, h);    % Channel path (convolution)
      yt = awgn(y_channel,EbN0,'measured');  
      y = yt(N_g+1:N_ofdm);       % Remove CP
      Y = fft(y);                           % FFT
k=1:N_p; Est_LS(k) = Y(pilot_loc(k))./X_p(k); % LS channel estimation
      Est_HLS = interpolate(Est_LS,pilot_loc,N_fft,'linear');
      h_estLS = ifft(Est_HLS); h_DFT = h_estLS(1:channel_length); 
      Est_HDFT = fft(h_DFT,N_fft); % DFT-based channel estimation
         MSE_LSi   =  MSE_LSi+ (H-Est_HLS )*(H-Est_HLS )';
         MSE_DFTi  =  MSE_DFTi+ (H-Est_HDFT)*(H-Est_HDFT)';
      end
      MSE_LS(i)=MSE_LSi;   MSE_DFT(i)=MSE_DFTi;
end
   MSE_LS = MSE_LS/(N_fft*N_sym);
   MSE_DFT = MSE_DFT/(N_fft*N_sym);
figure(1), semilogy(EbN0s',MSE_LS,'-x', EbN0s',MSE_DFT,'-d')
legend('LS-linear','LS-linearDFT')
xlabel ('Eb/N0 (dB)') ;
ylabel ( 'MSE')


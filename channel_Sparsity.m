M_H_range = [1 8 64];
M_V_range = [64 8 1];

%azimuth angle to be considered whe computing
varphiRangeDeg = (-90:0.1:90);
varphiRange = varphiRangeDeg*(pi/180);

% define channel properties
varphiDeg = [30 -20 40]; %azimuth angles of the treu channel path
varphi= varphiDeg*(pi/180);
theta = 0; % common elevation angle
%Normalized wavelength
lambda =1;

% now let's go through all three cases

for k = 1:length(M_H_range)
    M_H = M_H_range(k); %extracting number of antennas per row
    M_V = M_V_range(k); %extracting number of antennas per column
    d_H = 0.5*lambda;
    d_V = 0.5*lambda;
    
    %defining antenna geometry
    M = M_H*M_V; %total numberof antennas
    U = zeros(3,M); % creating matrix containing positions of antennas
    
    i =@(m)mod(m-1,M_H); %horizontal index
    j = @(m)floor((m-1)/M_H); % vertical index
    
    
    for m = 1:M
        U(:,m) = [0; i(m)*d_H; j(m)*d_V];
    end
    % Now compute the array response for various directions
    hRange = zeros(M,length(varphiRange));
    for n = 1:length(varphiRange)
        hRange(:,n) = functionSpartialSignature3DLoS(U,varphiRange(n),theta,lambda);
    end
    
    % channel computation by adding the array response for three different
    % path
    h = functionSpatialSignature3DLoS(U,varphi(1),theta,lambda);
    for n = 2:length(varphi)
        h = h + 0.5*functionSpartialSignature3DLoS(U,varphi(n),theta,lambda);
    end
    
    
    %compute how the received signal power is distributed over different
    %azimuth directions thi is channelgain that will be plotted
    
    gains = abs(h*hRange).^2;
    
    %plot the results
    
    figure(1);
    hold on; box on;
    plot(varphiRangeDeg, 10*log10(gains),'lineWidth',1);
    xlabel('Angleof arrival');
    ylabel('Signal gain [dB]');
end

legend('$1 \times 64$','$8 \times 8$', '$64 \times $1', 'Interpreter','latex','Location','southWest');

for n =1:length(varphiDeg)
    plot(varphiDeg(n)*[1 1],[-20 50],'k');
end
xlim([-90 90]);
ylim([-20 50]);
set(gca,'fontsize',14);
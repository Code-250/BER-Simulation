x =1;
Trial = 1000;
for EbN0 = 0:1:20
    linear_EbN0 = 10^(EbN0/10);
    nvar = 1/(linear_EbN0);
    err1 = 0;
    err2 = 0;
    err3 = 0;
        for trial = 1:Trial
            n1 = sqrt(nvar/2)*randn;
            n2= sqrt(nvar/2)*randn;
            h1 = sqrt(0.5)* abs(randn + j*randn);
            h2 = sqrt(0.5)* abs(randn + j*randn);
            
            y1 = x*h1+n1;
            y2 = x*h2*n2;
            y_equal = 0.5*(y1+y2);
            
            a1 = (abs(h1))^2;
            a2 = (abs(h2))^2;
            y_maximal = x*(a1*h1 + a2*h2)+a1*n1 + a2*n2;
            
            p1 = chi2rnd(4);
            p2 = chi2rnd(4);
            as1 = p1*(abs(h1))^2;
            as2 = p2*(abs(h2))^2;
            if as1 >= as2
                y_selection = x*(as1*h1)+as1*n1;
            end
            if as1 < as2
                y_selection = x*(as1 *h2)+as2*n2;
            end
            
            if y_equal <0
                err1 = err1 +1;
            end
            if y_maximal < 0
                err2 = err2 +1;
            end
            if y_selection < 0
                err3 = err2 + 1;
            end
        end
        BER1(EbN0+1) = err1/(Trial);
        BER2(EbN0+1) = err2/(Trial);
        BER3(EbN0+1) = err3/(Trial);
end
%plot simulations
figure
EbN0 = 0:1:20; %changed from 10
mu = 10.^(EbN0./10);
ber_theory = (1/2)*(1 - sqrt(mu ./ (mu + 1)));
semilogy(EbN0,BER1,'r*-',EbN0,BER2,'b--o',EbN0,BER3,'c-o',EbN0,ber_theory,'b'); %plot EG BER vs ebN0
legend('EG','NR','SC','theory');
xlabel('EbN0(dB)');
ylabel('bit error rate');
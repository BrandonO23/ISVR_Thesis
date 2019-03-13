clc
clear

load BafData.mat
load Meshiter.mat
% load matlab.mat
% load test3.mat


%%
deg = -175:5:180;
freq = D10(:,1);

iter = 1;
for i = 1:2:length(deg).*2
    D10C(:,iter) = complex(D10(:,i+1),D10(:,i+2));
    D20C(:,iter) = complex(D20(:,i+1),D20(:,i+2));
    D30C(:,iter) = complex(D30(:,i+1),D30(:,i+2));
    D40C(:,iter) = complex(D40(:,i+1),D40(:,i+2));
    D50C(:,iter) = complex(D50(:,i+1),D50(:,i+2));
    D60C(:,iter) = complex(D60(:,i+1),D60(:,i+2));
    
    iter = 1 + iter;
end

%%
fr = Abscissa;
iter = 1;
for i = 1:2:length(DD)*2
    M5000C(:,iter) = complex(M5000_60(:,i),M5000_60(:,1+i));
    M8000C(:,iter) = complex(M8000_60(:,i),M8000_60(:,1+i));
    M10000C(:,iter) = complex(M10000_60(:,i),M10000_60(:,1+i));
    M12000C(:,iter) = complex(M12000_60(:,i),M12000_60(:,1+i));
    M14000C(:,iter) = complex(M14000_60(:,i),M14000_60(:,1+i));
    M16000C(:,iter) = complex(M16000_60(:,i),M16000_60(:,1+i));
    M20000C(:,iter) = complex(M20000_60(:,i),M20000_60(:,1+i));
    iter = iter + 1;
end


%%
% 
% hold on
% plot(deg,20*log10((abs(D10C(end,:)))))
% plot(deg,20*log10((abs(D20C(end,:)))))
% plot(deg,20*log10((abs(D30C(end,:)))))
% plot(deg,20*log10((abs(D40C(end,:)))))
% plot(deg,20*log10((abs(D50C(end,:)))))
% plot(deg,20*log10((abs(D60C(end,:)))))
% 
% legend('10','20','30','40','50','60')
% grid on

%% Polar over freq
% for i = 1:size(D60C,1)
%     polar(deg.*pi./180,abs(D60C(i,:))./max(abs(D60C(i,:))))
%     title(num2str(freq(i)))
%     grid on
%     pause(.05)
% end

%%
% m1 = (max(abs(D60C(48,:))));
% rad = deg.*pi./180;
% 
% polar(deg.*pi./180,(abs(D60C(109,:))./max(abs(D60C(109,:)))))
% hold on
% polar(deg.*pi./180,abs(D60C(140,:))./max(abs(D60C(140,:))),'--')
% polar(deg.*pi./180,abs(D60C(155,:))./max(abs(D60C(155,:))),'-.')
% 
% legend('5K','8K','10K')

%% Calc DI

DI10 = calcDI(D10C,freq,DD);
DI20 = calcDI(D20C,freq,DD);
DI30 = calcDI(D30C,freq,DD);
DI40 = calcDI(D40C,freq,DD);
DI50 = calcDI(D50C,freq,DD);
DI60 = calcDI(D60C,freq,DD);

%%

DI8c = calcDI(M8000C,Abscissa,DD);
DI10c = calcDI(M10000C,Abscissa,DD);
DI12c = calcDI(M12000C,Abscissa,DD);
DI14c = calcDI(M14000C,Abscissa,DD);
DI16c = calcDI(M16000C,Abscissa,DD);
DI20c = calcDI(M20000C,Abscissa,DD);



%%
figure(1)
subplot(2,3,1)
plot(2.*pi.*.005.*freq./344,DI10);
subplot(2,3,2)
plot(2.*pi.*.01.*freq./344,DI20);
subplot(2,3,3)
plot(2.*pi.*.015.*freq./344,DI30);
subplot(2,3,4)
plot(2.*pi.*.02.*freq./344,DI40);
subplot(2,3,5)
plot(2.*pi.*.025.*freq./344,DI50);
subplot(2,3,6)
plot(2.*pi.*.03.*freq./344,DI60);

s = ['10';'20';'30';'40';'50';'60'];

for i = 1:6
    subplot(2,3,i)
    title([s(i,:) 'mm'])
    grid on
    ylabel('DI (dB)')
    ylim([-6 8])
    xlabel('ka')
end
    

%%
figure(2)
hold on
plot(2.*pi.*.03.*Abscissa./344,DI8c);
plot(2.*pi.*.03.*Abscissa./344,DI10c);
plot(2.*pi.*.03.*Abscissa./344,DI12c);
plot(2.*pi.*.03.*Abscissa./344,DI14c);
plot(2.*pi.*.03.*Abscissa./344,DI16c);
plot(2.*pi.*.03.*Abscissa./344,DI20c);

legend('8k','10k','12k','14k','16k','20k');

grid on
ylabel('DI (dB)')
title('mesh size')
ylim([-26 10])
xlabel('ka')
    



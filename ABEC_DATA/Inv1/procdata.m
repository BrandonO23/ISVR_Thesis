clc
clear
load inv1data
Freq = D10mm(:,1);


%%



for i = 1:length(degrees)
    D10C(:,i) = complex(D10mm(:,i+1),D10mm(:,i+2));
end
% ppol(D10C,degrees,Freq)

for i = 1:length(Freq)
    D = [D10C(i,1:15),D10C(i,19:end)]';
    B = [D10C(i,16:18)]';
    AC(i) = 10*log10((abs(B'*B).*21)./(abs(D'*D).*3));
end
semilogx(Freq,AC)



%%

for i = 1:length(degrees)
    D20C(:,i) = complex(D20mm(:,i+1),D20mm(:,i+2));
end
% ppol(D20C,degrees,Freq)


%%
for i = 1:length(degrees)
    D30C(:,i) = complex(D30mm(:,i+1),D30mm(:,i+2));
end
% ppol(D30C,degrees,Freq)

for i = 1:length(Freq)
    D = [D30C(i,1:11),D30C(i,13:end)]';
    B = D30C(i,12);
    AC(i) = 10*log10((abs(B'*B).*23)./(abs(D'*D)));
end


plot(AC)















%%
function ppol(DD,deg,f)
    for i = 1:2:201
        polarplot(deg.*pi/180,20*log10(abs(DD(i,:)./.00002)));
        title(['Freq = ' num2str(f(i))])
        rlim([min(20*log10(abs(DD(i,:)./.00002)))-2 max(20*log10(abs(DD(i,:)./.00002)))+2])
        pause(.5)
        
    end
end


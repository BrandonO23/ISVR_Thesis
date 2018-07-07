function DI = calcDI(x,freq,deg)
    Zeroind = find(deg == 0);
    M = abs(x(:,Zeroind));
    D = abs([x(:,1:Zeroind-1),x(:,Zeroind+1:end)]);
    for i = 1:length(freq)
        DI(i) = 10.*log10((M(i).^2./mean(D(i,:).^2)));
    end
end


%ACC_BEM

clc
clear
load D20
Freq = L20(:,1);

for i = 1:length(degr)
    L20C(:,i) = complex(L20(:,i+1),L20(:,i+2));
    R20C(:,i) = complex(R20(:,i+1),R20(:,i+2));
end

for f = 1:length(Freq)
    % Setup Acoustic variables
    omega = 2*pi*f;      % Angular frequency 
    c = 344;             % Speed of sound
    lambda = c./f;       % Wavelength
    rho = 1.225;         % Density of air
    k = (2*pi)./lambda;  % Wave number
    
    Zd = [[L20C(f,1:11),L20C(f,13:end)];[R20C(f,1:11),R20C(f,13:end)]].';
    Zb = [L20C(f,12),R20C(f,12)];
    nD = size(Zd,1);
    nB = 1;
    
    Rd = (Zd'*Zd);      
    Rb = (Zb'*Zb);
    
    %% Without regularization
    con(f) = cond(Rd);
    [V,D1] = eig(Rd\Rb);

    [md, idx] = max(diag(D1));
    q = V(:,idx);
    lam1 = sqrt(1./(q'*Rb*q));
    q = lam1.*q;
    
    Ref = 1j*omega*rho*exp(-1i*k.*1)./(4*pi*1); 
    qmono = mean(Zb*q)./Ref;
    
    AE(f) = 10*log10((q'*q)./((qmono'*qmono)));
    AC(f) = 10*log10((abs(q'*Rb*q.*nD))./(abs(q'*Rd*q.*nB)));
    
end
%%
figure(1)
% subplot(1,2,1)
semilogx(Freq,AC),title('Acoustic Contrast')
grid
% subplot(1,2,2)
% semilogx(Freq,AE),title('Array Effort')
% grid
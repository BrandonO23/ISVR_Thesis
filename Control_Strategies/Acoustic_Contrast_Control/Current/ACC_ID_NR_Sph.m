% Acoustic Contrast - Indirect Method - No Regularization - Sphere of
% microphones 

clc 
clear
freq = 10.^(2:.01:4);   % Frequency Vectors
rad = 1;                % Radius of Sphere
inc = 15;               % Degree Increments
fb = .005;              % Gives thickness of baffle/2
dd = .005;              % Distance from left to right of loudspeakers

%% Source positions in meters [x,y,z], can take any number of control sources 
Cs = [fb  dd 0;
      fb -dd 0;
     -fb  dd 0;
     -fb -dd 0];
  
[bpos,dpos,Sph,deg,nX,nY,nZ] = crossBD(rad,inc);

%% Freq Loop
for n = 1:length(freq)

    % Setup Acoustic variables
    omega = 2*pi*freq(n);      % Angular frequency 
    c = 344;                   % Speed of sound
    lambda = c./freq(n);       % Wavelength
    rho = 1.225;               % Density of air
    k = (2*pi)./lambda;        % Wave number

    % Length 
    for i = 1:size(Cs,1)
        D(:,i) = sqrt((dpos(:,1)-Cs(i,1)).^2 + (dpos(:,2)-Cs(i,2)).^2 + (dpos(:,3)-Cs(i,3)).^2);
        B(:,i) = sqrt((bpos(:,1)-Cs(i,1)).^2 + (bpos(:,2)-Cs(i,2)).^2 + (bpos(:,3)-Cs(i,3)).^2);
        Gd(:,i) = 1j*omega*rho*exp(-1i*k.*D(:,i))./(4*pi*D(:,i));
        Gb(:,i) = 1j*omega*rho*exp(-1i*k.*B(:,i))./(4*pi*B(:,i));
    end

    
    Lb = size(Gb,1);    % Number of control points in Bright zone
    Ld = size(Gd,1);    % Number of control points in Dark zone

    % Bright and Dark correletion matrices of acoustic trasnfer functions
    Rd = (Gd'*Gd)./Ld;      
    Rb = (Gb'*Gb)./Lb;
    
    Bright = 1;

    % Without regularization
    con(n) = cond(Rd);
    [V,D1] = eig(Rd\Rb);

    [md, idx] = max(diag(D1));
    q = V(:,idx);

    % Find lambda 
    lam1 = sqrt(Bright./(q'*Rb*q));
    q = lam1.*q;

    % Build single monopole reference for Array Effort
    Ref = 1j*omega*rho*exp(-1i*k.*rad)./(4*pi*rad); 
    qmono = mean(Gb*q)./Ref;
    
    AE(n) = 10*log10((q'*q)./((qmono'*qmono)/2));
    AC(n) = 10*log10((abs(q'*Rb*q))./(abs(q'*Rd*q)));
    
     
%   VISfield(Cs,q,f) 
%   pause(.1)

end
%%
figure(1)
subplot(1,2,1)
semilogx(freq,AC),title('Acoustic Contrast')
ylim([0 11])
grid
subplot(1,2,2)
semilogx(freq,AE),title('Array Effort')
ylim([-6 22])
grid

% figure(2)
% sphplot(Cs,dpos,bpos)

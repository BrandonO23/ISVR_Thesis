%Draft 1
clc 
clear
freq = 10.^(1:.005:4);
iter = 1;
rad = .5;

%% Source positions in meters [x,y], can take any number of control sources 
Cs = [ 0.04 0  0;
       0 0 0; 
      -0.04 0  0];
  
[Sph,deg,nX,nY,nZ] = evenSph(rad,15);
figure(1)
evenSph(rad,15);


for f = freq
% Setup Acoustic variables
omega = 2*pi*f;      % Angular frequency 
c = 344;             % Speed of sound
lambda = c./f;       % Wavelength
rho = 1.225;         % Density of air
k = (2*pi)./lambda;  % Wave number
  
%% Bright
Zind = find(nZ(:,1) == 0);
bind = find(deg == 0);
bpos = [nX(Zind,bind),0,0];
  
%% Dark
[~,indx]=ismember(bpos,Sph,'rows');
dpos = [Sph(1:indx(1)-1,:,:);Sph(indx(1)+1:end,:,:)];
  
%% Length 
for i = 1:size(Cs,1)
    D{i} = sqrt((dpos(:,1)-Cs(i,1)).^2 + (dpos(:,2)-Cs(i,2)).^2 + (dpos(:,3)-Cs(i,3)).^2);
    B{i} = sqrt((bpos(:,1)-Cs(i,1)).^2 + (bpos(:,2)-Cs(i,2)).^2 + (bpos(:,3)-Cs(i,3)).^2);
    Gd(:,i) = 1j*omega*rho*exp(-1i*k.*D{i})./(4*pi*D{i});
    Gb(:,i) = 1j*omega*rho*exp(-1i*k.*B{i})./(4*pi*B{i});
end

%% Build Full G matrix
G = [Gb;Gd];
a = [ones(size(Gb,1),1);zeros(size(Gd,1),1)];

%% Solve using PM
q = (G'*G)\G'*a;

%% Build single monopole reference for Array Effort
Ref = 1j*omega*rho*exp(-1i*k.*rad)./(4*pi*rad); 
qmono = mean(Gb*q)/Ref;

% Array effort and Acoustic Contrast
AE(iter) = 10*log10((q'*q)./((qmono'*qmono)));


%%
beta = 0;

while (AE(iter) >= 20) 
    q = (G'*G + beta*eye(3))\G'*a;
    beta = beta + 1;
    % Array effort and Acoustic Contrast
    AE(iter) = 10*log10((q'*q)./((qmono'*qmono)/2));
end

% Bright and Dark correletion matrices of acoustic transfer functions
Rd = (Gd'*Gd);      
Rb = (Gb'*Gb);
Lb = size(Gb,1);    % Number of control points in Bright zone
Ld = size(Gd,1);    % Number of control points in Dark zone

% Scale bright zone to 1 Pa 
% lam1 = sqrt(1./(q'*Rb*q));
% q = lam1.*q;

AC(iter) = 10*log10((Ld.*real(q'*Rb*q))./(Lb.*real(q'*Rd*q)));
% SPL(iter) = 20*log10(abs(mean(Gb*q))/.00002);

iter = iter + 1;
end
%%
figure(2)
subplot(1,2,1)
semilogx(freq,AC),title('Acoustic Contrast')
ylim([2 14])
grid
subplot(1,2,2)
semilogx(freq,AE),title('Array Effort')
ylim([-10 60])
grid

%%


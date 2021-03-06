%% Acoustic Contrast - Direct Method - No Regularization - No Array Effor Control
% Uses only one bright point on a sphere of sensors

clc 
clear
freq = 10.^(2:.0005:4);
iter = 1;
for f = freq
% Setup Acoustic variables
omega = 2*pi*f;      % Angular frequency 
c = 344;             % Speed of sound
lambda = c./f;       % Wavelength
rho = 1.225;         % Density of air
rad = .5;
k = (2*pi)./lambda;  % Wave number

%% Source positions in meters [x,y], can take any number of control sources 
Cs = [ 0.02 0  0;
      -0.02 0  0];
  
[Sph,deg,nX,nY,nZ] = evenSph(rad,15);
  
  
%% Bright
Zind = find(nZ(:,1) == 0);
bind = find(deg == 0);
bpos = [nX(Zind,bind),0,0];
  
%% Dark Points
[~,indx]=ismember(bpos,Sph,'rows');
dpos = [Sph(1:indx(1)-1,:,:);Sph(indx(1)+1:end,:,:)];

  
%% Length 
for i = 1:size(Cs,1)
    D(:,i) = sqrt((dpos(:,1)-Cs(i,1)).^2 + (dpos(:,2)-Cs(i,2)).^2 + (dpos(:,3)-Cs(i,3)).^2);
    B(:,i) = sqrt((bpos(:,1)-Cs(i,1)).^2 + (bpos(:,2)-Cs(i,2)).^2 + (bpos(:,3)-Cs(i,3)).^2);
    Gd(:,i) = 1j*omega*rho*exp(-1i*k.*D(:,i))./(4*pi*D(:,i));
    Gb(:,i) = 1j*omega*rho*exp(-1i*k.*B(:,i))./(4*pi*B(:,i));
end

%%
Db = 0.1;           % Constant real value for Dark zone
Lb = size(Gb,1);    % Number of control points in Bright zone
Ld = size(Gd,1);    % Number of control points in Dark zone

% Bright and Dark correletion matrices of acoustic trasnfer functions
Rd = (Gd'*Gd);      
Rb = (Gb'*Gb);

%% eigenvector & eigenvalues from Equation 3
[V,D1] = eig((Rd)\(Rb));

%% Max eigenvector corresponding to max eigenvalue
md = max(diag(D1));
N = find(diag(D1) == md);
q = V(:,N);

% Find lambda that satisfies Equation 4
lam1 = sqrt(Db./(q'*Rd*q));
q = lam1.*q;

% q = Rd\Gb';


%% Build single monopole reference for Array Effort
Ref = 1j*omega*rho*exp(-1i*k.*rad)./(4*pi*rad); 
qmono = mean(Gb*q)./Ref;


 
%% Array effort and Acoustic Contrast
AE(iter) = 10*log10((q'*q)./((qmono'*qmono)/2));
AC(iter) = 10*log10((real(q'*Rb*q.*Ld))./(real(q'*Rd*q.*Lb)));


iter = iter + 1;
end
%%
subplot(1,2,1)
semilogx(freq,AC),title('Acoustic Contrast')
% ylim([2 14])
grid
subplot(1,2,2)
semilogx(freq,AE),title('Array Effort')
% ylim([-10 60])
grid

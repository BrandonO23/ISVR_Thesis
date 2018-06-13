%% Acoustic Contrast - Indirect Method - No Regularization - 
% Uses only one bright point on a sphere of sensors

clc 
clear
freq = 10.^(2:.001:4);
iter = 1;
rad = 1;
%% Source positions in meters [x,y], can take any number of control sources 
Cs = [.09 0 0;
       .03 0 0;
      -.03 0 0;
      -.09 0 0];
  
[Sph,deg,nX,nY,nZ] = evenSph(rad,15);
  
  
%% Bright
Zind = find(nZ(:,1) == 0);
bpos(1,:) = [nX(Zind,find(deg == 0)), 0,0];
bpos(2,:) = [nX(Zind,find(deg == 15)) ,nY(Zind,find(deg == 15)),0];
bpos(3,:) = [nX(Zind,find(deg == 360-15)),nY(Zind,find(deg == 360-15)),0];
bpos(4,:) = [nX(5,find(deg == 0)),0,nZ(5,find(deg == 0))];
bpos(5,:) = [nX(7,find(deg == 0)),0,nZ(7,find(deg == 0))];


%% Dark Points
[~,indx]=ismember(bpos,Sph,'rows');
indx = sort(indx);
dpos = [Sph(1:indx(1)-1,:,:);
        Sph(indx(1)+1:indx(2)-1,:,:);
        Sph(indx(2)+1:indx(3)-1,:,:);
        Sph(indx(3)+1:indx(4)-1,:,:);
        Sph(indx(4)+1:indx(5)-1,:,:);
        Sph(indx(5)+1:end,:,:)];

%% Freq Loop
for f = freq

    % Setup Acoustic variables
    omega = 2*pi*f;      % Angular frequency 
    c = 344;             % Speed of sound
    lambda = c./f;       % Wavelength
    rho = 1.225;         % Density of air
    k = (2*pi)./lambda;  % Wave number

    %% Length 
    for i = 1:size(Cs,1)
        D(:,i) = sqrt((dpos(:,1)-Cs(i,1)).^2 + (dpos(:,2)-Cs(i,2)).^2 + (dpos(:,3)-Cs(i,3)).^2);
        B(:,i) = sqrt((bpos(:,1)-Cs(i,1)).^2 + (bpos(:,2)-Cs(i,2)).^2 + (bpos(:,3)-Cs(i,3)).^2);
        Gd(:,i) = 1j*omega*rho*exp(-1i*k.*D(:,i))./(4*pi*D(:,i));
        Gb(:,i) = 1j*omega*rho*exp(-1i*k.*B(:,i))./(4*pi*B(:,i));
    end

    %%
    Lb = size(Gb,1);    % Number of control points in Bright zone
    Ld = size(Gd,1);    % Number of control points in Dark zone

    % Bright and Dark correletion matrices of acoustic trasnfer functions
    Rd = (Gd'*Gd);      
    Rb = (Gb'*Gb);

    %% Without regularization
    con(iter) = cond(Rd);
    [V,D1] = eig(Rd\Rb);

    [md, idx] = max(diag(D1));
    q = V(:,idx);


    %% Build single monopole reference for Array Effort
    Ref = 1j*omega*rho*exp(-1i*k.*rad)./(4*pi*rad); 
    qmono = mean(Gb*q)./Ref;
    
    AE(iter) = 10*log10((q'*q)./((qmono'*qmono)/2));
    AC(iter) = 10*log10((real(q'*Rb*q.*Ld))./(real(q'*Rd*q.*Lb)));
    
    iter;
    iter = iter + 1;
end
%%
figure(1)
subplot(1,2,1)
semilogx(freq,AC),title('Acoustic Contrast')
ylim([0 12])
grid
subplot(1,2,2)
semilogx(freq,AE),title('Array Effort')
% ylim([-5 40])
grid

figure(2)
scatter3(dpos(:,1),dpos(:,2),dpos(:,3),'k','fill'),hold on
scatter3(bpos(:,1),bpos(:,2),bpos(:,3),'k')
scatter3(Cs(:,1),Cs(:,2),Cs(:,3),150,'+','r')

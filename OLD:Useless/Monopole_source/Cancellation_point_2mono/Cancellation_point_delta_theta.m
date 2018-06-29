%% Cancellation at point in far field

clc
clear

f = 1000;             % Single frequency
omega = 2*pi*f;      % Angular frequency 
c = 344;             % Speed of sound
lambda = c./f;       % Wavelength
rho = 1.225;         % Density of air
k = 2*pi./lambda;    % Wave number
qp = .0005;               % Volume Velocity


% To make a bit easier the field is limited to first quadrant
delta = .01;
rx = -1:delta:1;       % Radius in x
ry = -1:delta:1;        % Radius in y
lx = length(rx);     % length of rx
ly = length(ry);     % length of ry
d = .05;              % Distance from origin to monopole source
rcan = 1;           % Distance to cancellation


[X, Y] = meshgrid(rx,ry);        % Meshgrid from rx and ry
Z = sqrt((X).^2 + (Y).^2);       % Field from Origin
Zp = sqrt((X+d).^2 + (Y).^2);    % Field seen by Source 1
Zc = sqrt((X-d).^2 + (Y).^2);    % Field seen by Source 2


for theta = 1:5:180
Zcan = ceil([rcan*cos(theta*pi/180),rcan*sin(theta*pi/180)]./delta);
[zr,zc] = find(Z < 1e-5);
zc = zc + Zcan(1);
zr = zr + 44; 
distp = Zp(zr,zc);
distc = Zc(zr,zc);


a = (omega*rho./(4*pi*distc))^2;
b = (omega.*rho./(4*pi)).^2.*(1/(distp.*distc)).*qp.*exp((-1i*k).*(distp-distc));
c = (omega*rho./(4*pi*distp))^2;
qc = -inv(a)*b;


pp = 1j*omega*rho*qp*exp(-1i*k.*Zp)./(4*pi*Zp); 
pc = 1j*omega*rho*qc*exp(-1i*k.*Zc)./(4*pi*Zc); 
p = pp + pc;

surf(rx,ry,(abs(p)),'edgecolor', 'none')
caxis([-1 1])
hold on
scatter3(0-d,0,2,'o','linewidth',2,'MarkerFaceColor','k','MarkerEdgeColor','k')
scatter3(0+d,0,2,'o','linewidth',2,'MarkerFaceColor','k','MarkerEdgeColor','k')
hold off

colormap('jet')
view(0,90)
colorbar
xlabel('Meters'),ylabel('Meters')
title(['Cancellation at ', num2str(rcan),' meters, \theta = ',num2str(theta),'.   Frequency = ',num2str(f),' Hz'])
% M(theta) = getframe;
pause(.001)
end

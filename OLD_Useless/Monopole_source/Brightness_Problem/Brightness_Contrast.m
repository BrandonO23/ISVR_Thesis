%% Brightness Problem
clc
clear
it = 1;


f = 5000;            % Single frequency
omega = 2*pi*f;      % Angular frequency 
c = 344;             % Speed of sound
lambda = c./f;       % Wavelength
rho = 1.225;         % Density of air
k = 2*pi./lambda;    % Wave number
Jo = (5e-5)^2;      % square modulus volume Velocity Contraint
ref = .00002;

delta = .01;
rx = -1:delta:1;                 % x
ry = -1:delta:1;                  % y
[X, Y] = meshgrid(rx,ry);        % Meshgrid from rx and ry


% Bright Zone
b = [.5,0];
len = .1/delta;
wid = .1/delta;
row = find(abs(10000000000*(rx-b(1)))<1);
col = find(abs(10000000000*(ry-b(2)))<1);

len = row:row-1 + .1/delta;
wid = col:col-1 + .1/delta;
tS = length(len)*length(wid);
[b1, b2] = meshgrid(wid,len);
bZ1 = reshape(b1,tS,1);
bZ2 = reshape(b2,tS,1);

% Dark Zone
d = [-.5,.2];
len = .1/delta;
wid = .1/delta;
row = find(abs(10000000000*(rx-d(1)))<1);
col = find(abs(10000000000*(ry-d(2)))<1);

len = row:row-1 + .1/delta;
wid = col:col-1 + .1/delta;
tS = length(len)*length(wid);
[d1, d2] = meshgrid(wid,len);
dZ1 = reshape(d1,tS,1);
dZ2 = reshape(d2,tS,1);


% Source positions in meters
Cs = [.005 0;
        0,0;
     -.005 0];
  
% Pre Allocation
l = size(Cs,1);                  
MeshZ = cell(l,1);
green = cell(l,1);
Gb = zeros(tS,l);
Gd = zeros(tS,l);
for i = 1:l
    MeshZ{i} = sqrt((X-Cs(i,1)).^2 + (Y-Cs(i,2)).^2);
    green{i} = 1j*omega*rho*exp(-1i*k.*MeshZ{i})./(4*pi*MeshZ{i});
    for j = 1:tS
        Gb(j,i) = green{i}(bZ1(j),bZ2(j));
        Gd(j,i) = green{i}(dZ1(j),dZ2(j));
    end
end
Rd = sum(Gb)'*sum(Gb)./tS;
Rb = sum(Gd)'*sum(Gd)./tS;

[V,D] = eig(inv(Rb)*Rd);
qc = V(:,2);

p = zeros(length(ry),length(rx));
lam = Jo/(qc'*qc);
q = lam.*qc;

for i = 1:l
    p = p + green{i}.*q(i);
end


lam = Jo/(qc'*qc);
q = lam.*qc;

surf(rx,ry,(abs(p)),'edgecolor', 'none')
colormap('jet')
caxis([-.5 .5])
view(0,90)
co = colorbar;
co.Label.String = 'SPL (dB ref @ 20 \muPa)';
xlabel('Meters'),ylabel('Meters')
hold on
scatter3(b(1),b(2),10000,'o','linewidth',2,'MarkerFaceColor','k','MarkerEdgeColor','k')
hold off

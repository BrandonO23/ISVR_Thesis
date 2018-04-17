%Draft 1
clc 
clear

f = 100;
omega = 2*pi*f;      % Angular frequency 
c = 344;             % Speed of sound
lambda = c./f;       % Wavelength
rho = 1.225;         % Density of air
k = 2*pi./lambda;    % Wave number
Jo = (.8*5e-5)^2;      % square modulus volume Velocity Contraint

delta = .01;
rx = -1:delta:1;                 % x
ry = 0:delta:2;                  % y
[X, Y] = meshgrid(rx,ry);        % Meshgrid from rx and ry
radius = 1;
Meshtest = sqrt((X).^2 + (Y).^2);
lag = 2;

%%
% Measurement Points
pos = zeros(180/5,2);
int = 1;
deg = 0:5:180;
for i = deg

    xd = radius*cos(i*pi/180);
    yd = radius*sin(i*pi/180);

    dx  =  find(round(X(1,:),2) == round(xd,2));
    dy  =  find(round(Y(:,1),2) == round(yd,2));
    
    pos(int,:) = [dx,dy]; int = int + 1;
end
%%
% Bright 90 degrees
for bdeg = 90:90
% bdeg = 45;
bpos = pos(deg == bdeg,:);

% Dark everywhere else
dpos = [pos(1:find(deg == bdeg)-1,:); pos(find(deg == bdeg)+1:end,:)];

%%
% Source positions in meters
Cs = [.0094 0; 
      .0047 0;
     -.0047 0;
     -.0094 0]; 
  
l = size(Cs,1);                  
MeshZ = cell(l,1);
green = cell(l,1);
p = zeros(length(ry),length(rx));


for i = 1:l
    MeshZ{i} = sqrt((X-Cs(i,1)).^2 + (Y-Cs(i,2)).^2);
    green{i} = 1j*omega*rho*exp(-1i*k.*MeshZ{i})./(4*pi*MeshZ{i});
    for j = 1:length(dpos)
        Gd(j,i) = green{i}(dpos(j,2),dpos(j,1));
    end
    Gb(i) = green{i}(bpos(2),bpos(1));
end

[V,D] = eig((Gd'*Gd + lag*eye(4))\(Gb'*Gb));
q = V(:,find(diag(D) == max(diag(D))));


 for i = 1:l
    p = p + green{i}.*q(i);
 end
%%

surf(rx,ry,real(p),'edgecolor', 'none')
colormap('jet')
caxis([0 5])
view(0,90)
colorbar
xlabel('Meters'),ylabel('Meters')
hold on
for i = 1:length(dpos)
    scatter3(X(1,dpos(i,1)),Y(dpos(i,2),1),10000,'o','linewidth',2,'MarkerFaceColor','w','MarkerEdgeColor','w')
end
scatter3(X(1,bpos(1,1)),Y(bpos(1,2),1),10000,'x','linewidth',2,'MarkerFaceColor','w','MarkerEdgeColor','w')
hold off


pause(.5)
end
%%
% % Define Bright Zone
% by =  1; % meter in y
% bx = -0.2000; % meters wide/2
% zr  =   find(round(X(1,:),2) == 0);
% xindL = find(round(X(1,:),2) == bx);
% xindR = find(round(X(1,:),2) == -bx);
% yind = find(round(Y,2) == by);
% stopind = r(1) - yind(1); 
% 
% % Bright Vector
% BV = zeros(stopind+1,1);
% for i = xindL:xindR
%     holder = (yind(i):(yind(i) + stopind)).';
%     BV = [BV; holder];
% end
%         




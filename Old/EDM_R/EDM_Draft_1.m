%Draft 1
clc 
clear
it = 1;
f = 100:10:500;
qit = zeros(length(f),1);

%%
iter = 1;
% for f = 1000:1000:10000
f = 1000;
omega = 2*pi*f;      % Angular frequency 
c = 344;             % Speed of sound
lambda = c./f;       % Wavelength
rho = 1.225;         % Density of air
k = 2*pi./lambda;    % Wave number

delta = .01;
rx = -1:delta:1;                 % x
ry = -1:delta:1;                  % y
[X, Y] = meshgrid(rx,ry);        % Meshgrid from rx and ry
radius = .5;
Meshtest = sqrt((X).^2 + (Y).^2);


%%
% Measurement Points

int = 1;
deg = 0:15:360 - 15;
pos = zeros(length(deg),2);
for i = deg

    xd = radius*cos(i*pi/180);
    yd = radius*sin(i*pi/180);

    dx  =  find(round(X(1,:),2) == round(xd,2));
    dy  =  find(round(Y(:,1),2) == round(yd,2));
    
    pos(int,:) = [dx,dy]; int = int + 1;
end
%%
% Bright 90 degrees
% for d = 0:5:160
% d = 35;
bdeg = 75:15:105;
bind = find(deg == bdeg(1)):find(deg == bdeg(end));
bpos = pos(bind,:);

% Dark everywhere else
if bind(1) == 1
    dpos = pos(bind(end)+1:end,:);
else
    dpos = [pos(1:bind(1)-1,:); pos(bind(end)+1:end,:)];
end
%%
% Source positions in meters
Cs = [.004 0;
      .03 0;
     -.03 0;
     -.004 0;];
  
l = size(Cs,1);                  
MeshZ = cell(l,1);
green = cell(l,1);
p = zeros(length(ry),length(rx));

Gb = zeros(size(bpos,1),l);
Gd = zeros(size(dpos,1),l);

for i = 1:l
    MeshZ{i} = sqrt((X-Cs(i,1)).^2 + (Y-Cs(i,2)).^2);
    green{i} = 1j*omega*rho*exp(-1i*k.*MeshZ{i})./(4*pi*MeshZ{i});
    for j = 1:length(dpos)
        Gd(j,i) = green{i}(dpos(j,2),dpos(j,1));
    end
    for j = 1:size(bpos,1)
        Gb(j,i) = green{i}(bpos(j,2),bpos(j,1));
    end
end

%%

E = (.8*5e-5)^2; 
Rd = (Gd'*Gd)./size(Gd,1);
Rb = (Gb'*Gb)./size(Gd,1);
beta = 10^7;
alpha = 1;
% for alpha = a
[V,D] = eig((Gb'*Gb - alpha.*Gd'*Gd));

md = max(diag(D));
q = V(:,find(diag(D) == md));
q1 = md*V(:,find(diag(D) == md));


% while q'*q >= E
% 
%     [V,D] = eig((Gd'*Gd + beta*eye(size(Cs,1)))\(Gb'*Gb));
%     md = max(diag(D));
%     q = md*V(:,find(diag(D) == md));
% %     test(iter) = q'*q;
% %     iter = iter + 1;
%     beta = beta + 100000000;
% 
% end
AE(iter) = q1'*q1;
Pre(iter) = 20*log10(abs(Gb*q1)/.00002);
iter = iter + 1;

% plot(a,Pre)
 %%
 
 for i = 1:l
    p = p + green{i}.*q(i);
 end

%%

surf(rx,ry,abs(p),'edgecolor', 'none')
colormap('jet')
caxis([0 10])
view(0,90)
colorbar
xlabel('Meters'),ylabel('Meters')
title(['F = ' num2str(f)]);
hold on
for i = 1:length(dpos)
    scatter3(X(1,dpos(i,1)),Y(dpos(i,2),1),10000,'o','linewidth',2,'MarkerFaceColor','w','MarkerEdgeColor','w')
end
for i = 1:size(bpos,1)
        scatter3(X(1,bpos(i,1)),Y(bpos(i,2),1),10000,'x','linewidth',2,'MarkerFaceColor','w','MarkerEdgeColor','w')
end
hold off
pause(.5)
% end

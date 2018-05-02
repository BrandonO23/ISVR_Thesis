%Draft 1
clc 
clear
freq = 10.^(1:.005:4);
iter = 1;
for f = freq
% Setup Acoustic variables
omega = 2*pi*f;      % Angular frequency 
c = 344;             % Speed of sound
lambda = c./f;       % Wavelength
rho = 1.225;         % Density of air
k = (2*pi)./lambda;  % Wave number

% Build Mesh Grid
delta = .01;                       % Descrete Spacing 
rx = -.5:delta:.5;                   % x 
ry = -.5:delta:.5;                   % y
rz = -.5:delta:.5;                   % z
[X, Y, Z] = meshgrid(rx,ry, rz);     % X Y and Z to create meshgrid
radius = .3;                       % Radius of circular control points

MeshRef = sqrt((X).^2 + (Y).^2 + (Z).^2);   % Reference Mesh distance


%% Measurement Points

intM = 1;                       % iter value for pos array
deg = 0:15:360 - 15;            % array of degree points on circle
pos = zeros(length(deg),2);     % Preallocate position vector 

% Loop that finds the position of each point on the circular array and puts
% it in the pos vector
for i = deg

    xd = radius*cos(i*pi/180);  % x component
    yd = radius*sin(i*pi/180);  % y component

    % Finds position on Meshgrid using X and Y
    dx  =  find(round(X(1,:),2) == round(xd,2)); 
    dy  =  find(round(Y(:,1),2) == round(yd,2));
    
    % Inserts that position in pos vector
    pos(intM,:) = [dx,dy]; 
    intM = intM + 1;
end

%% Bright Points
bdeg = 0:0;     % Vector that gives the position on circular array

% Finds indices of bright points, and pulls them from pos
bind = find(deg == bdeg(1)):find(deg == bdeg(end));
bpos = pos(bind,:);

%% Dark points
% Puts the rest of the points as dark points
if bind(1) == 1
    dpos = pos(bind(end)+1:end,:);
else
    dpos = [pos(1:bind(1)-1,:); pos(bind(end)+1:end,:)];
end

%% Source positions in meters [x,y], can take any number of control sources 
Cs = [ 0.02 0;
      -0.02 0];
  
% Preallocate l meshgrids, greens functions, and pressure matrices/vectors 
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

%% Build Full G matrix
G = [Gb;Gd];
a = [ones(size(Gb,1),1);zeros(size(Gd,1),1)];


%% Solve using PM
q = (G'*G)\G'*a;



%% Solve for total field, using superposition
for i = 1:l
   p = p + green{i}.*q(i);
end

%% Build single monopole refeernece for Array Effort
Gr = zeros(size(bpos,1),1);
for j = 1:size(bpos,1)
     greenR = 1j*omega*rho*exp(-1i*k.*MeshRef)./(4*pi*MeshRef);
     Gr(j) = greenR(bpos(j,2),bpos(j,1));
end
avgG = mean(Gr); 
qmono = mean(Gb*q)/avgG;

% Bright and Dark correletion matrices of acoustic trasnfer functions
Rd = (Gd'*Gd);      
Rb = (Gb'*Gb);
Lb = size(Gb,1);    % Number of control points in Bright zone
Ld = size(Gd,1);    % Number of control points in Dark zone
 
% Array effort and Acoustic Contrast
AE(iter) = 10*log10((q'*q)./((qmono'*qmono)));
AC(iter) = 10*log10((Ld.*real(q'*Rb*q))./(Lb.*real(q'*Rd*q)));


%%
% figg = figure(1);
% surf(rx,ry,20*log10(abs(p)./.00002),'edgecolor', 'none')
% colormap('jet')
% caxis([40 100])
% view(0,90)
% colorbar
% xlabel('Meters'),ylabel('Meters')
% title(['F = ' num2str(f)]);
% hold on
% for i = 1:length(dpos)
%     scatter3(X(1,dpos(i,1)),Y(dpos(i,2),1),10000,'o','linewidth',2,'MarkerFaceColor','w','MarkerEdgeColor','w')
% end
% for i = 1:size(bpos,1)
%         scatter3(X(1,bpos(i,1)),Y(bpos(i,2),1),10000,'x','linewidth',2,'MarkerFaceColor','w','MarkerEdgeColor','w')
% end
% hold off
% pause(.01)
 iter = iter + 1;
end
%%
subplot(1,2,1)
semilogx(freq,AC),title('Acoustic Contrast')
grid
subplot(1,2,2)
semilogx(freq,AE),title('Array Effort')
grid

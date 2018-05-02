% Build Mesh Grid
delta = .01;                       % Descrete Spacing 
rx = -10:delta:10;                   % x 
ry = -10:delta:10;                   % y
[X, Y] = meshgrid(rx,ry);          % X and Y to create meshgrid
MeshRef = sqrt((X).^2 + (Y).^2);   % Reference Mesh distance


omega = 2*pi*[100,200];      % Angular frequency 
c = 344;             % Speed of sound
lambda = c./[100,200];       % Wavelength
rho = 1.225;         % Density of air
k = (2*pi)./lambda;  % Wave number


% Source positions in meters [x,y], can take any number of control sources 
Cs = [ 5 0;
      -5 0];
  
q = [.0015;
     .0007];
% Preallocate l meshgrids, greens functions, and pressure matrices/vectors 
l = size(Cs,1);                                 
MeshZ = cell(l,1);
green = cell(l,1);
p = zeros(length(ry),length(rx));

for i = 1:l
    MeshZ{i} = sqrt((X-Cs(i,1)).^2 + (Y-Cs(i,2)).^2);
    green{i} = 1j*omega(i)*rho*exp(-1i*k(i).*MeshZ{i})./(4*pi*MeshZ{i});
end

for i = 1:l
   p = p + green{i}.*q(i);
end

figg = figure(1);
surf(rx,ry,real(p),'edgecolor', 'none')
colormap('jet')
caxis([-.1 .1])
view(0,90)
colorbar
xlabel('Meters'),ylabel('Meters')

title('Pressure field Real(Pa), Left = 200 Hz, Right = 100 Hz')

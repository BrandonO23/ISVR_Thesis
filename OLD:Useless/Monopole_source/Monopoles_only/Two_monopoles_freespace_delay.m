%% Two monopoles changine delay
clc
clear

f = 3000;            % Single frequency
omega = 2*pi*f;      % Angular frequency 
c = 344;             % Speed of sound
lambda = c./f;       % Wavelength
rho = 1.225;         % Density of air
k = 2*pi./lambda;    % Wave number
q = .05;             % Volume Velocity
rx = -1:.01:1;       % Radius in x
ry = -1:.01:1;       % Radius in y
lx = length(rx);     % length of rx
ly = length(ry);     % length of ry
d = .1;
delay = 0.001:.01:1;

[X, Y] = meshgrid(rx,ry);        % Meshgrid from rx and ry
Z1 = sqrt((X).^2 + (Y-d).^2);    % Distance to point from origin
Z2 = sqrt((X).^2 + (Y+d).^2);    % Distance to point from origin


for i = 1:length(delay)
    pz1 = exp(-1i*k.*Z1);               % First Pressure field 
    pz2 = exp(-1i*k.*Z2+delay(i));      % Second Pressure fieldPressure
    pz = pz1 + pz2;                     % If linear
    p = pz;
    surf(rx,ry,(real(p)),'edgecolor', 'none')
    colormap('jet')
    view(0,90)
    colorbar
    title(['delay = ', num2str(delay(i))])
    pause(.001)
    xlabel('Meters'),ylabel('Meters')    
end
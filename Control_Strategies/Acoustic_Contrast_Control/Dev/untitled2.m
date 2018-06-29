% Acoustic Contrast - Indirect Method - No Regularization - Circle of
% microphones 

clc 
clear
freq = 10.^(2:.01:4);   % Frequency Vectors

%% Measurement Points

intM = 1;                       % iter value for pos array
deg = 0:5:360 - 5;            % array of degree points on circle
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
bdeg = 0;     % Vector that gives the position on circular array

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


%% Build Mesh Grid
delta = .01;                       % Descrete Spacing 
rx = -1:delta:1;                   % x 
ry = -1:delta:1;                   % y
[X, Y] = meshgrid(rx,ry);          % X and Y to create meshgrid
radius = 1;                        % Radius of circular control points
MeshRef = sqrt((X).^2 + (Y).^2);   % Reference Mesh distance

for n = 1:length(freq)

    % Setup Acoustic variables
    omega = 2*pi*freq(n);      % Angular frequency 
    c = 344;                   % Speed of sound
    lambda = c./freq(n);       % Wavelength
    rho = 1.225;               % Density of air
    k = (2*pi)./lambda;        % Wave number

    % Length 
    for i = 1:size(Cs,1)
        D(:,i) = sqrt((dpos(:,1)-Cs(i,1)).^2 + (dpos(:,2)-Cs(i,2)).^2 + (dpos(:,3)-Cs(i,3)).^2);
        B(:,i) = sqrt((bpos(:,1)-Cs(i,1)).^2 + (bpos(:,2)-Cs(i,2)).^2 + (bpos(:,3)-Cs(i,3)).^2);
        Gd(:,i) = 1j*omega*rho*exp(-1i*k.*D(:,i))./(4*pi*D(:,i));
        Gb(:,i) = 1j*omega*rho*exp(-1i*k.*B(:,i))./(4*pi*B(:,i));
    end

    
    Lb = size(Gb,1);    % Number of control points in Bright zone
    Ld = size(Gd,1);    % Number of control points in Dark zone

    % Bright and Dark correletion matrices of acoustic trasnfer functions
    Rd = (Gd'*Gd)./Ld;      
    Rb = (Gb'*Gb)./Lb;
    
    Bright = 1;

    % Without regularization
    con(n) = cond(Rd);
    [V,D1] = eig(Rd\Rb);

    [md, idx] = max(diag(D1));
    q = V(:,idx);

    % Find lambda 
    lam1 = sqrt(Bright./(q'*Rb*q));
    q = lam1.*q;

    % Build single monopole reference for Array Effort
    Ref = 1j*omega*rho*exp(-1i*k.*rad)./(4*pi*rad); 
    qmono = mean(Gb*q)./Ref;
    
    AE(n) = 10*log10((q'*q)./((qmono'*qmono)/2));
    AC(n) = 10*log10((abs(q'*Rb*q))./(abs(q'*Rd*q)));
    
     
%   VISfield(Cs,q,f) 
%   pause(.1)

end


%Draft 1
clc 
clear
freq = 10.^(3:.05:4.3);
iter = 1;
itit = 1;
% for dgg = 0:5:360-1
for f = freq
    % Setup Acoustic variables

    omega = 2*pi*f;      % Angular frequency 
    c = 344;             % Speed of sound
    lambda = c./f;       % Wavelength
    rho = 1.225;         % Density of air
    k = (2*pi)./lambda;  % Wave number

    % Build Mesh Grid
    delta = .01;                      % Descrete Spacing 
    rx = -1:delta:1;                 % x 
    ry = -1:delta:1;                 % y
    [X, Y] = meshgrid(rx,ry);          % X and Y to create meshgrid
    radius = 1;                        % Radius of circular control points
    MeshRef = sqrt((X).^2 + (Y).^2);   % Reference Mesh distance


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

    %% From Reference Document

    Db = 0.1;           % Constant real value for Dark zone
    Lb = size(Gb,1);    % Number of control points in Bright zone
    Ld = size(Gd,1);    % Number of control points in Dark zone

    % Bright and Dark correletion matrices of acoustic trasnfer functions
    Rd = (Gd'*Gd);      
    Rb = (Gb'*Gb);

    % eigenvector & eigenvalues from Equation 3
    [V,D] = eig((Rd)\(Rb));

    % Max eigenvector corresponding to max eigenvalue
    md = max(diag(D));
    q = V(:,find(diag(D) == md));

    % Find lambda that satisfies Equation 4
    lam1 = sqrt(Db./(q'*Rd*q));
    q = lam1.*q;

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

    %% Array effort and Acoustic Contrast
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
%     scatter3(X(1,dpos(i,1)),Y(dpos(i,2),1),10000,'o','linewidth',2,'MarkerFaceColor','k','MarkerEdgeColor','k')
% end
% for i = 1:size(bpos,1)
%         scatter3(X(1,bpos(i,1)),Y(bpos(i,2),1),10000,'x','linewidth',2,'MarkerFaceColor','k','MarkerEdgeColor','k')
% end
% hold off
% pause(.01)
iter = iter + 1;
end
% figure(1)
% subplot(1,2,1)
% semilogx(freq,(AC)),title('Acoustic Contrast')
% grid
% subplot(1,2,2)
% semilogx(freq,AE),title('Array Effort')
% grid
% pause(.01)
% avAC(itit) = mean(AC);
% avAE(itit) = mean(AE);
% itit = itit + 1;
% AC = [];
% AE = [];
% iter = 1;
% end
%%
subplot(1,2,1)
semilogx(freq,(AC)),title('Acoustic Contrast')
ylim([0 5])
grid
subplot(1,2,2)
semilogx(freq,AE),title('Array Effort')
grid

%%
% figure(2)
% subplot(1,2,1)
% plot(0:5:360-1,avAC)
% ylabel('Average Acoustic Contrast'),xlabel('Degrees')
% axis tight,grid
% title('2 sources - 4 cm apart')
% subplot(1,2,2)
% plot(0:5:360-1,avAE)
% ylabel('Average Array Effort'),xlabel('Degrees')
% axis tight,grid
% title('2 sources - 4 cm apart')
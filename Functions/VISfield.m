function VISfield(Cs,q,f)

    % Build Mesh Grid
    delta = .01;                       % Descrete Spacing 
    rx = -2:delta:2;                   % x 
    ry = -2:delta:2;                   % y
    [X, Y] = meshgrid(rx,ry);          % X and Y to create meshgrid


    omega = 2*pi*f;      % Angular frequency 
    c = 344;             % Speed of sound
    lambda = c./f;       % Wavelength
    rho = 1.225;         % Density of air
    k = (2*pi)./lambda;  % Wave number

    % Preallocate l meshgrids, greens functions, and pressure matrices/vectors 
    l = size(Cs,1);                                 
    MeshZ = cell(l,1);
    green = cell(l,1);
    p = zeros(length(ry),length(rx));

    for i = 1:l
        MeshZ{i} = sqrt((X-Cs(i,1)).^2 + (Y-Cs(i,2)).^2);
        green{i} = 1j*omega*rho*exp(-1i*k.*MeshZ{i})./(4*pi*MeshZ{i});
    end

    for i = 1:l
       p = p + green{i}.*q(i);
    end

    surf(rx,ry,abs(p),'edgecolor', 'none')
    colormap('jet')
    caxis([-1 1])
    view(0,90)
    colorbar
    xlabel('Meters'),ylabel('Meters')
    title(['Pressure field Real(Pa), F = ', num2str(f) ])

end

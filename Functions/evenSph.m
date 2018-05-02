function [Sph,deg,nX,nY,nZ] = evenSph(Radius,degrees)
% evenSph builds a sphere with points evenly spaced by degrees and with
% Radius
    
    % Degree Vector
    deg = 0:degrees:360 - degrees;
    
    % Build full sphere with, 24 for 0:15:360-1;
    [X,Y,Z] = sphere(size(deg,2));

    % Get rid of the center points, and last row
    Xtemp = X(2:end-1,1:end-1);
    Ytemp = Y(2:end-1,1:end-1);
    Ztemp = Z(2:end-1,1:end-1);

    % Dummy Var for Loop
    it = 1;

    % loop that only takes even numbers
    for i = 2:2:size(X,1)-2
        nX(it,:)= Radius.*Xtemp(i,:);
        nY(it,:)= Radius.*Ytemp(i,:);
        nZ(it,:)= Radius.*Ztemp(i,:);
        it = it +1;
    end

    Sph = [nX(:),nY(:),nZ(:);0,0,Radius;0,0,-Radius];

    if nargout == 0
        scatter3(Sph(:,1),Sph(:,2),Sph(:,3))
    end
end
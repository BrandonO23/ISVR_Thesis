function [bpos,dpos,Sph,deg,nX,nY,nZ] = crossBD(rad,inc)


% Builds Sphere of of equally spaced points with radius rad and degree 
% increment 15 
[Sph,deg,nX,nY,nZ] = evenSph(rad,inc);
  
  
%% Bright Zone points
%  Set as a cross at the 0 degree point

Zind = find(nZ(:,1) == 0);
bpos(1,:) = [nX(Zind,find(deg == 0)), 0,0];
bpos(2,:) = [nX(Zind,find(deg == 15)) ,nY(Zind,find(deg == 15)),0];
bpos(3,:) = [nX(Zind,find(deg == 360-15)),nY(Zind,find(deg == 360-15)),0];
bpos(4,:) = [nX(5,find(deg == 0)),0,nZ(5,find(deg == 0))];
bpos(5,:) = [nX(7,find(deg == 0)),0,nZ(7,find(deg == 0))];


%% Dark Zone Points

[~,indx]=ismember(bpos,Sph,'rows');
indx = sort(indx);
dpos = [Sph(1:indx(1)-1,:,:);
        Sph(indx(1)+1:indx(2)-1,:,:);
        Sph(indx(2)+1:indx(3)-1,:,:);
        Sph(indx(3)+1:indx(4)-1,:,:);
        Sph(indx(4)+1:indx(5)-1,:,:);
        Sph(indx(5)+1:end,:,:)];
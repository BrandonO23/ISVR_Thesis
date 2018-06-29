function [bpos,dpos,Sph,deg,nX,nY,nZ] = circleBD(rad,inc)


% Builds Sphere of of equally spaced points with radius rad and degree 
% increment 15 
[Sph,deg,nX,nY,nZ] = evenSph(rad,inc);
  
  
%% Bright Zone points
%  Set as a cross at the 0 degree point

Nd = length(deq)/2;

for i = 1:3
%     bpos(i,:) =  
end
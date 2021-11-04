function PlotEarth
%-------------------------------------------------------------------------%
%
% plotEarth.m plots the Earth inside a 3d graphic.
%
%-------------------------------------------------------------------------%
% CALLED FUNCTIONS:  
%  (none)
%
% ------------------------------------------------------------------------%
% CONTRIBUTORS:
%  Gian Marco Paldino
%  Gabriele Palumbo
%  Matteo Zeni
%
%-------------------------------------------------------------------------%
% VERSIONS: 
%  26/12/2020: first version
%
%-------------------------------------------------------------------------%

C = imread('map.jpg');
theta = 0;
[x,y,z] = ellipsoid(0,0,0, 6378.135, 6378.135,6356.750,1E2);
surf(x,y,z,circshift(flip(C),[0,ceil(size(C,2)/360*theta)]),'FaceColor','texturemap','EdgeColor','none');
axis equal;
end

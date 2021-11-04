function [thvect,rr,vv] = plotOrbit(a,e,i,OM,om,thi,thf,dth,mu,color,width)
%-------------------------------------------------------------------------%
%
% plotOrbit.m propagates and plots the unperturbed orbit starting from 
% Keplerian coordinates.
% 
%-------------------------------------------------------------------------%
% PROTOTYPE:
%  [thvect,rr,vv] = plotOrbit(a,e,i,OM,om,thi,thf,dth,mu,color,width)
%
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  a            [1]    Semi-major axis                      [km]
%  e            [1]    Eccentricity                         [-]
%  i            [1]    Inclination                          [rad]
%  OM           [1]    RAAN                                 [rad]
%  om           [1]    Pericenter anomaly                   [rad]
%  thi          [1]    Initial true anomaly                 [rad]
%  thf          [1]    Final true anomaly                   [rad]
%  dth          [1]    True anomaly step                    [rad]
%  mu           [1]    Standard gravitational parameter     [km^3/s^2]
%  color        [char] Plot color                           [-]
%  width        [char] Plot linewidth                       [-]
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  thvect       [1xlength(thvect)] True anomaly vector      [rad]
%  rr           [3xlength(thvect)] Position matrix          [km]
%  vv           [3xlength(thvect)] Velocity matrix          [km/s]
%  Actual plot
%
%-------------------------------------------------------------------------%
% CALLED FUNCTIONS:  
%  kep2car.m
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

% Theta vector, position matrix, velocity matrix are initializated

thvect = (thi:dth:thf);

rr = zeros(3,length(thvect));
vv = zeros(3,length(thvect));

% Position vector and velocity vector are calculated for each theta

for n = 1:length(thvect)
    
    [rrs,vvs] = kep2car(a,e,i,OM,om,thvect(n),mu);
    rr(:,n) = rrs;
    vv(:,n) = vvs;
end

% Positions are plotted into 3D graphic 

plot3(rr(1,:),rr(2,:),rr(3,:),color,'linewidth',width);
grid on
hold on

end
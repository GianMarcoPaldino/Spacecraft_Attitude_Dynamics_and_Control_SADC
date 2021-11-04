function [rr, vv] = kep2car(a,e,i,OM,om,th,mu)
%-------------------------------------------------------------------------%
%
% kep2car.m trasforms Keplerian parameters to cartesian coordinates.
%
%-------------------------------------------------------------------------%
% PROTOTYPE:
%  [rr, vv] = kep2car(a, e, i, OM, om, th, mu)
% 
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS:
%  a             [1]     Semi-major axis                      [km]
%  e             [1]     Eccentricity                         [-]
%  i             [1]     Inclination                          [rad]
%  OM            [1]     RAAN                                 [rad]
%  om            [1]     Argument of periapsis                [rad]
%  th            [1]     True anomaly                         [rad]
%  mu            [1]     Standard gravitational parameter     [km^3/s^2]
%                        of the primary body
%
%-------------------------------------------------------------------------%
% OUTPUT ARGUMENTS:
%  rr            [3x1]   Position vector                      [km]
%  vv            [3x1]   Velocity vector                      [km/s]
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
%  26-12-2020: First version
%
%-------------------------------------------------------------------------%

    % Semi-latus rectum
    p = a*(1-e^2); 
    
    % Position vector norm
    r = p/(1+e*cos(th));
        
    % Position and velocity vectors in Perifocal Coordinates [3x1]
    rr_pf = r.*[cos(th); sin(th); 0];
    
    vv_pf = (sqrt(mu/p))*[-sin(th); e+cos(th); 0];
    
    if i == 0 
        warning('The orbit has 0 inclination, thus there ascending node is not defined. The angle of apse line from X axis is OM + om')
    end
    
    % First rotation matrix (around k versor)
    R_OM = [cos(OM) , sin(OM) , 0;...
           -sin(OM) , cos(OM) , 0;...
           0        , 0       , 1];

    % Second rotation matrix (around i' versor)
    R_i  = [1      , 0       ,     0 ;...
            0      , cos(i)  , sin(i);...
            0      , -sin(i) , cos(i)];   
        
    % Third rotation matrix (around k'' versor)
    R_om = [cos(om) , sin(om) ,     0 ;...
            -sin(om), cos(om) ,     0 ;...
            0       , 0       ,     1 ];
    
    % Rotation matrix from ECI to PF
    T_eci2pf = R_om * R_i * R_OM;
    
    % Rotation matrix from PF to ECI
    T_pf2eci = T_eci2pf.';
   
    % Position and velocity vectors rotated in ECI 
    rr = T_pf2eci*rr_pf;
    vv = T_pf2eci*vv_pf;
       
end
   
    
    


%-------------------------------------------------------------------------%
%                                                                         %                             
%                       SPACECRAFT ATTITUDE DYNAMICS                      %
%                              A.Y. 2020/2021                             %
%                                                                         %
%                        GIAN MARCO PALDINO - 968731                      % 
%                                                                         %
%-------------------------------------------------------------------------%
clc
clear all 
addpath(genpath(pwd))
%-------------------------------------------------------------------------%
% Orbit Keplerian parameters 
a = 7000;                    %  Semi-major axis        [km]
e = 0.0100;                  %  Eccentricity           [-]
i = deg2rad(20);             %  Inclination            [rad]
OM = deg2rad(55);            %  RAAN                   [rad]
om = deg2rad(70);            %  Argument of perigee    [rad]
theta = deg2rad(0);          %  True anomaly           [rad]

% Orbit's initial Keplerian elements vector 
kep = [a; e; i; OM; om; theta]; 

% Earth's standard gravitational parameter [km^3/s^2]
mu = 398600;

% Orbital period [s]
T_orb = 2*pi*sqrt(a^3/mu);

% Inertia moments [kg * m^2]
Ix = 0.135;
Iy = 0.140;
Iz = 0.252;

I = [Ix Iy Iz];

% Initial angular velocity 
w_0 = deg2rad([8 13 4]);

% Initial attitude matrix
A_BN_0 = eye(3);

% Load IGRF g and h quasi-normalized coefficients 
load('IGRF_coefficients.mat')

%% S/C geometrical properties

% Dimensions [m]
% Main body
w = 0.2263;
l = 0.366;
h = 0.1; 
% Solar panels 
w_SP = 0.325;
h_SP = 0.223;

% Normals, Surface Coefficients, Areas:
N_body = 6; % # main body surfaces
N_SP = 4; % # solar array surfaces

N=zeros(N_body+N_SP,3); % [10 x 3] Matrix containing (row-wise) the normal of each surface in Body Frame
A=zeros(N_body+N_SP,1); % [10 x 1] Column vector Vector containing the area of each surface
d=zeros(N_body+N_SP,3); % [10 x 3] Matrix containing (row-wise) the position vector of each surface wrt to the center 

n1 = [1 0 0];  n2 = [0 1 0];  n3 = [0 0 1]; % outward normal unit vectors

% Main Body 
N(1,:) = n1;    N(2,:) = n2;     N(3,:) = n3;
N(4,:) = -n1;   N(5,:) = -n2;    N(6,:) = -n3;
c_s_body = .5;     
c_d_body = .1;
A(1:N_body) = [w*h l*h w*l w*h l*h w*l]; % [m^2]

% Solar Panels
N(7,:) = n1;    N(8,:) = -n1;   N(9,:) = n1;    N(10,:) = -n1;
c_s_SP = .8;    c_d_SP = .1;
A(N_body+1:N_body+N_SP) = [w_SP*h_SP w_SP*h_SP w_SP*h_SP w_SP*h_SP]; % [m^2]

% Reference Frame centered @ Geometric Center of the body:
% Coordinates of each surface's center from the Geometric Center [m]
d(1,:) = [l/2 0 0];  d(2,:) = [0 w/2 0];  d(3,:) = [0 0 h/2];
d(4,:) = [-l/2 0 0]; d(5,:) = [0 -w/2 0]; d(6,:) = [0 0 -h/2];
d(7,:) = [l/2 (w/2 + w_SP/2) 0];   d(8,:)  = [l/2 (w/2 + w_SP/2) 0];
d(9,:) = [l/2 -(w/2 + w_SP/2) 0];   d(10,:) = [l/2 -(w/2 + w_SP/2) 0];

% Displacement of the Center of Gravity wrt Geometric Center [m]
d_cg = [0.017 0 0]'; 

%% Start Simulink 
open nocontrol

%% Perform simulation
sim nocontrol

%% Plot orbit
figure(1)
[thvect,rr,vv] = plotOrbit(a,e,i,OM,om,0,2*pi,0.001,398600,'r',2);
PlotEarth
xlabel('$r_x$ [km]','Interpreter','latex')
ylabel('$r_y$ [km]','Interpreter','latex')
zlabel('$r_z$ [km]','Interpreter','latex')

%% Define simulation time vector
t = out.tout;

%% Angular velocity plot
figure(2)
plot(t,rad2deg(out.w.Data),'linewidth',1)
legend('\omega_x','\omega_y','\omega_z')
grid on
xlim([0,2*T_orb])
xlabel('Time [s]')
ylabel('Angular velocity [deg/s]')

%% Environmental disturbance torques
figure(3)
subplot(2,2,1)
plot(t,out.T_GG)
grid on
legend('$T_{x}^{GG}$', '$T_y^{GG}$', '$T_{z}^{GG}$','interpreter','latex')
xlim([0,2*T_orb])
xlabel('Time [s]')
ylabel('Gravity Gradient torque [Nm]')

subplot(2,2,2)
plot(t,out.T_SRP)
grid on
legend('$T_{x}^{SRP}$', '$T_{y}^{SRP}$', '$T_{z}^{SRP}$','interpreter','latex')
xlim([0,2*T_orb])
xlabel('Time [s]')
ylabel('SRP torque [Nm]')

subplot(2,2,3)
plot(t,out.T_magn)
grid on
legend('$T_{x}^{magn}$', '$T_{y}^{magn}$', '$T_{z}^{magn}$','interpreter','latex')
xlim([0,2*T_orb])
xlabel('Time [s]')
ylabel('Magnetic torque [Nm]')

subplot(2,2,4)
plot(t,out.T_drag)
grid on
legend('$T_{x}^{drag}$', '$T_{y}^{drag}$', '$T_{z}^{drag}$','interpreter','latex')
xlim([0,2*T_orb])
xlabel('Time [s]')
ylabel('Atmospheric drag torque [Nm]')

figure(4)
plot(t,out.T_mag)
grid on
xlim([0,2*T_orb])
xlabel('Time [s]')
ylabel('Overall disturbance torque magnitude [Nm]')

%% Order of magnitude of disturbance torques

T_mag_ = [out.T_GG_mag out.T_SRP_mag out.T_magn_mag out.T_drag_mag];
T_ordersOfMagnitude = mean(abs(T_mag_));

h = figure;
set(h, 'Units', 'normalized', 'OuterPosition', [0.25 0.12 0.5 0.8]);
bar(T_ordersOfMagnitude)
set(gca,'YScale','log')
grid on
name = {'GG';'SRP';'Magnetic';'Drag'};
set(gca,'xticklabel',name)
ylabel('Average torque magnitude {[Nm]}')

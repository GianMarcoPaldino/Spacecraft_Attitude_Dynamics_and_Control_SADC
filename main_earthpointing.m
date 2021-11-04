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
w_0 = deg2rad([0 0 0]);

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
open Earth_pointing

%% Perform simulation
sim Earth_pointing

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
subplot(2,1,1)
plot(t,rad2deg(out.w.Data),t,rad2deg(out.w_est(:,:)),'--','linewidth',1)
legend('$\omega_x^{true}$','$\omega_y^{true}$','$\omega_z^{true}$','$\omega_x^{est}$','$\omega_y^{est}$','$\omega_z^{est}$','interpreter','latex')
grid on
xlim([0,T_orb])
xlabel('Time [s]')
ylabel('Angular velocity [deg/s]')

subplot(2,1,2)
plot(t,rad2deg(out.w_err(:,:)),'linewidth',1)
legend('$\omega_x^{est} - \omega_x^{true}$','$\omega_y^{est} - \omega_y^{true}$','$\omega_z^{est} - \omega_z^{true}$','interpreter','latex')
grid on
xlim([0,T_orb])
xlabel('Time [s]')
ylabel('Estimation error [deg/s]')


%% Angular velocity estimation error (average value)
w_err_mean = mean(abs(rad2deg(out.w_err(:,:))));

%% Attitude matrix
A_BN_11(1,:) = out.A_BN(1,1,:);
A_BN_12(1,:) = out.A_BN(1,2,:);
A_BN_13(1,:) = out.A_BN(1,3,:);
A_BN_21(1,:) = out.A_BN(2,1,:);
A_BN_22(1,:) = out.A_BN(2,2,:);
A_BN_23(1,:) = out.A_BN(2,3,:);
A_BN_31(1,:) = out.A_BN(3,1,:);
A_BN_32(1,:) = out.A_BN(3,2,:);
A_BN_33(1,:) = out.A_BN(3,3,:);

figure(3)
subplot(3,1,1)
plot(t,A_BN_11, t,A_BN_12, t,A_BN_13, t,A_BN_21, t,A_BN_22, t,A_BN_23, ...
    t,A_BN_31, t,A_BN_32,'b',t,A_BN_33,'r','linewidth',0.8)
grid on
xlim([0,T_orb])
xlabel('Time [s]')
ylabel('True DCM')
legend('$A_{B/N}^{11}$','$A_{B/N}^{12}$','$A_{B/N}^{13}$',...
    '$A_{B/N}^{21}$','$A_{B/N}^{22}$','$A_{B/N}^{23}$','$A_{B/N}^{31}$',...
    '$A_{B/N}^{32}$','$A_{B/N}^{33}$','interpreter','latex')

%Estimated attitude matrix and related estimation error
A_BN_svd_11(1,:) = out.A_BN_svd(1,1,:);
A_BN_svd_12(1,:) = out.A_BN_svd(1,2,:);
A_BN_svd_13(1,:) = out.A_BN_svd(1,3,:);
A_BN_svd_21(1,:) = out.A_BN_svd(2,1,:);
A_BN_svd_22(1,:) = out.A_BN_svd(2,2,:);
A_BN_svd_23(1,:) = out.A_BN_svd(2,3,:);
A_BN_svd_31(1,:) = out.A_BN_svd(3,1,:);
A_BN_svd_32(1,:) = out.A_BN_svd(3,2,:);
A_BN_svd_33(1,:) = out.A_BN_svd(3,3,:);


% Define attitude error matrix 
A_BN_err_11 = A_BN_svd_11 - A_BN_11;
A_BN_err_12 = A_BN_svd_12 - A_BN_12;
A_BN_err_13 = A_BN_svd_13 - A_BN_13;
A_BN_err_21 = A_BN_svd_21 - A_BN_21;
A_BN_err_22 = A_BN_svd_22 - A_BN_22;
A_BN_err_23 = A_BN_svd_23 - A_BN_23;
A_BN_err_31 = A_BN_svd_31 - A_BN_31;
A_BN_err_32 = A_BN_svd_32 - A_BN_32;
A_BN_err_33 = A_BN_svd_33 - A_BN_33;

subplot(3,1,2)
plot(t,A_BN_svd_11, t,A_BN_svd_12, t,A_BN_svd_13, t,A_BN_svd_21, ...
    t,A_BN_svd_22, t,A_BN_svd_23, t,A_BN_svd_31, t,A_BN_svd_32,'b',...
    t,A_BN_svd_33,'r')
legend('$\hat{A}_{B/N}^{11}$','$\hat{A}_{B/N}^{12}$',...
    '$\hat{A}_{B/N}^{13}$','$\hat{A}_{B/N}^{21}$','$\hat{A}_{B/N}^{22}$',...
    '$\hat{A}_{B/N}^{23}$','$\hat{A}_{B/N}^{31}$','$\hat{A}_{B/N}^{32}$',...
    '$\hat{A}_{B/N}^{33}$','interpreter','latex')
grid on
xlim([0,T_orb])
xlabel('Time [s]')
ylabel('Estimated DCM')

subplot(3,1,3)
plot(t,A_BN_err_11, t,A_BN_err_12, t,A_BN_err_13, t,A_BN_err_21, ...
    t,A_BN_err_22, t,A_BN_err_23, t,A_BN_err_31, t,A_BN_err_32,'b', ...
    t,A_BN_err_33, 'r')
legend('$\hat{A}_{B/N}^{11}-A_{B/N}^{11}$','$\hat{A}_{B/N}^{12}-A_{B/N}^{12}$',...
    '$\hat{A}_{B/N}^{13}-A_{B/N}^{13}$','$\hat{A}_{B/N}^{21}-A_{B/N}^{21}$',...
    '$\hat{A}_{B/N}^{22}-A_{B/N}^{22}$',...
    '$\hat{A}_{B/N}^{23}-A_{B/N}^{23}$','$\hat{A}_{B/N}^{31}-A_{B/N}^{31}$',...
    '$\hat{A}_{B/N}^{32}-A_{B/N}^{32}$',...
    '$\hat{A}_{B/N}^{33}-A_{B/N}^{33}$','interpreter','latex')
grid on
xlim([0,T_orb])
xlabel('Time [s]')
ylabel('DCM estimation error')

%% Desired attitude matrix and attitude error
% Desired attitude
A_LN_11(1,:) = out.A_LN(1,1,:);
A_LN_12(1,:) = out.A_LN(1,2,:);
A_LN_13(1,:) = out.A_LN(1,3,:);
A_LN_21(1,:) = out.A_LN(2,1,:);
A_LN_22(1,:) = out.A_LN(2,2,:);
A_LN_23(1,:) = out.A_LN(2,3,:);
A_LN_31(1,:) = out.A_LN(3,1,:);
A_LN_32(1,:) = out.A_LN(3,2,:);
A_LN_33(1,:) = out.A_LN(3,3,:);

figure(4)
subplot(2,1,1)
plot(t,A_LN_11, t,A_LN_12, t,A_LN_13, t,A_LN_21, t,A_LN_22, t,A_LN_23, ...
    t,A_LN_31, t,A_LN_32,'b',t,A_LN_33,'r','linewidth',0.8)
grid on
xlim([0,T_orb])
ylim([-1.1,1.1])
xlabel('Time [s]')
ylabel('Desired attitude')
legend('$A_{d}^{11}$','$A_{d}^{12}$','$A_{d}^{13}$',...
    '$A_{d}^{21}$','$A_{d}^{22}$','$A_{d}^{23}$','$A_{d}^{31}$',...
    '$A_{d}^{32}$','$A_{d}^{33}$','interpreter','latex')


% Attitude error
A_BL_11(1,:) = out.A_BL(1,1,:);
A_BL_12(1,:) = out.A_BL(1,2,:);
A_BL_13(1,:) = out.A_BL(1,3,:);
A_BL_21(1,:) = out.A_BL(2,1,:);
A_BL_22(1,:) = out.A_BL(2,2,:);
A_BL_23(1,:) = out.A_BL(2,3,:);
A_BL_31(1,:) = out.A_BL(3,1,:);
A_BL_32(1,:) = out.A_BL(3,2,:);
A_BL_33(1,:) = out.A_BL(3,3,:);

subplot(2,1,2)
plot(t,A_BL_11, t,A_BL_12, t,A_BL_13, t,A_BL_21, t,A_BL_22, t,A_BL_23, ...
    t,A_BL_31, t,A_BL_32,'b',t,A_BL_33,'c','linewidth',0.8)
grid on
xlim([0,T_orb])
ylim([-1.1,1.1])
xlabel('Time [s]')
ylabel('Desired attitude')
legend('$A_{e}^{11}$','$A_{e}^{12}$','$A_{e}^{13}$',...
    '$A_{e}^{21}$','$A_{e}^{22}$','$A_{e}^{23}$','$A_{e}^{31}$',...
    '$A_{e}^{32}$','$A_{e}^{33}$','interpreter','latex')

%% Attitude determination algorithms comparison

% Define error vectors
err_svd=out.err_svd;
err_alg = out.err_alg;

% Compute average values
err_svd_mean = mean(err_svd);
err_alg_mean = mean(err_alg);

% Plot results
figure(5)
plot(t,err_alg,t,err_svd)
legend('TRIAD','SVD')
xlim([10,T_orb])
grid on
xlabel('Time [s]')
ylabel('Attitude estimation error')

%% Overall control torque

u_id = zeros(length(t),3);
u_id_0 =zeros(length(t)-length(out.u_id),3);
u_id(1:length(u_id_0),:) = u_id_0;
u_id(length(u_id_0)+1:end,:) = out.u_id;
u_id(u_id==0) = nan;

figure(6)
u_cmg = out.u;
u_cmg(u_cmg==0) = nan;
subplot(2,1,2)
plot(t,u_cmg,'linewidth',0.8)
grid on
xlabel('Time [s]')
ylabel('Actual control torque [Nm]')
legend('$u_x$','$u_y$','$u_z$','interpreter','latex')

subplot(2,1,1)
plot(t,u_id)
grid on
xlabel('Time [s]')
ylabel('Ideal control law [Nm]')
legend('$u_{id}^x$','$u_{id}^y$','$u_{id}^z$','interpreter', 'latex')
magnifyOnFigure

%% Misalignment angle

mis_angle_x(1,:) = out.mis_angle_id(:,1,:);
mis_angle_y(1,:) = out.mis_angle_id(:,2,:);
mis_angle_z(1,:) = out.mis_angle_id(:,3,:);

figure(7)
subplot(3,1,1)
plot(t,mis_angle_x)
xlim([200,T_orb])
grid on
ylabel('X axis')

subplot(3,1,2)
plot(t,mis_angle_y)
xlim([200,T_orb])
grid on
ylabel('Y axis')

subplot(3,1,3)
plot(t,mis_angle_z)
xlim([200,T_orb])
grid on
xlabel('Time [s]')
ylabel('Z axis')

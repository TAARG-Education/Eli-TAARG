%% BladeSection_AngleOfAttack.m Test Case on UH-60

clc; close all; clear all;

%% Data
R = 8.17;            % Rotor radius [m/s]
c = 0.527;           % Chord [m]
A = pi*R^2;          % Disk area [m^2]
Cl_alpha = 5.73;     % Blade lift slope
Cd_mean = 0.01;      % Mean drag coefficiente
N = 4;               % Number of blades
sigma = 0.082;       % Blade solidity
omega_tip = 220.98;  % Tip speed [m/s]
omega = omega_tip/R; % Angular velocity [rad/s]
theta_tw = -18;      % Blade twist [deg]
W = 9979.0321;       % Weigth [kg]
W = 9.81*W;          % Weight [N]
h = 0;               % Altitude [m]             
Lock = 8.1936;       % Lock number
f = 3.376;           % Equivalent drag area [m^2]
X = 0;               % Rate of climb [m/s]
V_inf = 50;          % Asymptotic velocity [m]

r_segn = linspace(0.1925,1,100); 
psi = linspace(0,2*pi,100);
alpha_stall_up = 17;  % [deg]
alpha_stall_lo = -10; % [deg]

%% Solving articulated rotor
% Function needs to be modified for UH-60
[Tc,Hc,Yc,Qc,Pc,alfa,lambda,theta_0,beta_0,beta_1c,beta_1s] = ...
    Ndim_Coeff_Articulated_Rotor(V_inf,h,Lock,f,X);

%% Calculation
mu = V_inf*cos(convang(alfa,'deg','rad'))/(omega_tip);
beta = beta_0 + beta_1c*cos(psi) + beta_1s*sin(psi);
dbeta = -beta_1c*sin(psi) + beta_1s*cos(psi);
theta_tw = convang(-18,'deg','rad');
theta = theta_0 + theta_tw*r_segn;

[alpha_eff,Meff] = BladeSection_AngleOfAttack(r_segn,psi,mu,alfa,lambda, ...
    beta,dbeta,theta,h,omega_tip,alpha_stall_up,alpha_stall_lo,'on');

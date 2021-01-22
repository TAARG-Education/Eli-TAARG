%%% Test case
%%% Agusta Bell AB 206 Helicopter in forward flight
% Advance ratio = 0.12
% Descent angle = 20°


% User input
%%
clc;clear all;close all
global M R N c gamma theta_tw f rho      

M = 1120;                             % [Kg], mass
R = 5.1;                              % [m], rotorblade radius
N = 2;                                % number of blades
c = 0.34;                             % [m], chord length
gamma = 9;                            % Lock's number
theta_tw = -0.2313;                   % [rad], pitch gradient along the blade (linear variation is assumed)
f = 0.007;                            % [m^2] equivalent wet area
rho = 1.225;                          % [Kg/m^3], density (SML)

%% Function output
[T_C H_C alpha_deg omega V_inf]=autorot_performaces(0.12,20);
clc;clear all;close all
%% Test case

%% Agusta Bell AB 206 Helicopter in forward flight
% Advance ratio = 0.12
% Descent angle = 20°


% User input:

%Helicopter
M = 1120;                             % [Kg], mass
R = 5.1;                              % [m], rotorblade radius
N = 2;                                % number of blades
c = 0.34;                             % [m], chord length
gamma = 9;                            % Lock's number
theta_tw = -0.2313;                   % [rad], pitch gradient along the blade (linear variation is assumed)
f = 0.007;                            % [m^2] equivalent wet area

%Advance ratio
mu = 0.12;

%Descend angle
X = 20;

%% Function output
[T_C H_C alpha_deg omega V_inf]=autorot_performaces(mu,X,R,N,c,gamma,theta_tw,f,M);




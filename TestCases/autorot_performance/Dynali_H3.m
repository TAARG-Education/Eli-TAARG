clc;clear all;close all
%% Test case

%% Dynali H3 Sport Helicopter in forward flight

% User input:

%Helicopter
W = 500;                             % [Kg], mass
R = 7.14/2;                          % [m], rotorblade radius
N = 2;                               % number of blades
c = 0.18;                            % [m], chord length
gamma = 9;                           % Lock's number
theta_tw = -7.48;                    % [°], blade twisting (linear variation is assumed)
f = 0.007;                           % [m^2] equivalent wet area

%Advance ratio 
mu = 0.15;

%Descend angle
X = 20;

%% Function output
[T_C H_C alpha omega V_inf Theta_0]=autorot_performance(W,R,N,c,gamma,theta_tw,f,mu,X);

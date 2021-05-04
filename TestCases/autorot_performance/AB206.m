clc;clear all;close all
%Test case
%% Agusta Bell AB 206 Helicopter in forward flight
%  Advance ratio = 0.15
%  Descent angle = 20°


%User input:

%Helicopter
W = 1120;                             % [Kg], mass
R = 5.1;                              % [m], rotorblade radius
N = 2;                                % number of blades
c = 0.34;                             % [m], chord length
gamma = 9;                            % Lock's number
theta_tw = -13.2;                     % [°],blade twisting (linear variation is assumed)
f = 0.007;                            % [m^2] equivalent wet area

%Advance ratio
mu = 0.15;

%Descend angle
X = 25;

%% Function output
[T_C H_C alpha omega V_inf Theta_0]=autorot_performance(W,R,N,c,gamma,theta_tw,f,mu,X);




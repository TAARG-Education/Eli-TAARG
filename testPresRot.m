clc; clear; close all

Rt = 5.42;
Rr = Rt * .132;
thetaT = 13.5;
deltaTheta = 5;
omega = 240;
chord = .31;
N = 4;
Vc = 0;
rho = 1.225;
airfoil = 4412;


PresRot( Rr, Rt, thetaT, deltaTheta, omega, chord, N, Vc, rho, airfoil );
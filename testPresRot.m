clc; clear; close all

Rt = 2.286 / 2;
Rh = .178;
thetaT = 5;
deltaTheta = 0;
omega = 237;
chord = .38;
N = 2;
Vc = 0;
rho = 1.225;
airfoil = 0012;


[ Tc , Qc ] = ...
    PresRot( Rh, Rt, thetaT, deltaTheta, omega, chord, N, Vc, rho, airfoil )
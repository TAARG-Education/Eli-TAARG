clc; clear; close all;

chord = .192;
Rt = 3.5;
Rh = chord;
thetaHd = 10;
deltaThetad = 0;
omega = 45.72;

N = 4;
Vc = 0;
rho = 1.225;
airfoil = 0012;

[ Tc , Qc ] = ...
    RotPerf( Rh, Rt, thetaHd, deltaThetad, omega, chord, N, Vc, rho, airfoil )

Ts = Tc / .07
Qs = Qc / .07

%T = Tc * rho * omega^2 * pi * Rt^4;
%Q = Qc * rho * omega^2 * pi * Rt^5;

%Ct = T / ( rho * (omega/(2*pi))^2 * (2*Rt)^4 )
%Cq = Q / ( rho * (omega/(2*pi))^2 * (2*Rt)^5 )

%Cts = Ct / .07
%Cqs = Cq / .07
% TEST CASE 1

chord = .1905;
Rt = 6 * chord;
Rh = chord;
thetaHd = 5;
deltaThetad = 0;
rpm = 1250;
omega = rpm / 60 * 2 * pi;

N = 2;
Vc = 0;
rho = 1.225;
airfoil = 0012;

[ Tc , Qc ] = ...
    RotPerf( Rh, Rt, thetaHd, deltaThetad, omega, chord, N, Vc, rho, airfoil )

T = Tc * rho * omega^2 * pi * Rt^4;
Ct = T / ( rho * (rpm/60)^2 * (2*Rt)^4 )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST CASE 2
Rt = 11.9 / 2;
Rh = Rt * .2 ;
thetaHd = 9;
deltaThetad = 4;
omega = 240 / Rt;
chord = .385;
N = 4;
Vc = 0;
rho = 1.225;
airfoil = 0012;
%}

[ Tc , Qc ] = ...
    RotPerf( Rh, Rt, thetaHd, deltaThetad, omega, chord, N, Vc, rho, airfoil )

T = Tc * rho * omega^2 * pi * Rt^4;
Ct = T / ( rho * (rpm/59)^2 * (2*Rt)^4 )
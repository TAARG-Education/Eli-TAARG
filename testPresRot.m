clc; clear; close all;

Rt = convlength( 26.75, 'in', 'm' );
AR = 18.2;
chord = Rt / AR;
Rh = Rt * .148;
thetaHd = 8;
deltaThetad = 0;
omega = 45.72;
N = 4;
Vc = 0;
rho = 1.225;
airfoil = 0012;
sigma = chord * N / ( pi * Rt );

[ Tc , Qc ] = ...
    RotPerf( Rh, Rt, thetaHd, deltaThetad, omega, chord, N, Vc, rho, airfoil )

Ts = Tc / sigma
Qs = Qc / sigma

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
Ct = T / ( rho * (rpm/60)^2 * (2*Rt)^4 )
%% \main.m
%  \brief: this script tests the Ndim_Coeff_Articulated_Rotor.m function
%  \author: Marco Artiano, Luca Angelino
%  \version: 1.01
%
% Eli-TAARG is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public
% License as published by the Free Software Foundation; either
% version 3 of the License, or (at your option) any later version.
%
% Eli-TAARG is developed by the TAARG Educational organization for
% educational purposes only.
% Theoretical and Applied Aerodynamic Research Group - University of Naples Federico II.
%
% Eli-TAARG GitHub link: <https://github.com/TAARG-Education/Eli-TAARG>
%
% Eli-TAARG is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
% General Public License for more details.
% <http://www.gnu.org/licenses/>.
%
% ===============================================================================================
% |Name        : Main.m                                                                         |
% |Author      : Marco Artiano & Luca Angelino                                                  |
% |              University of Naples Federico II.                                              |
% |Version     : 1.0.4                                                                          |
% |Date        : 25/11/20                                                                       |
% |Modified    : 03/03/23 : Vincenzo Castrignano                                                |
% |Description : Determination of non-dimensional aerodynamics coefficients                     |
% |Reference   : 'Lezioni di AERODINAMICA DELL'ALA ROTANTE a.a. 2019-2020 - Renato Tognaccini'  |
% |Input       : Velocity V_inf [m/s], Altitude h [m], Lock number L [~],                       |            
% |              Equivalent drag area f [m^2], Rate of climb X [deg],                           |                                                             
% |              Rotor solidity sigma [~], Weight W [Kg], Blade radius R [m],                   |
% |              Rotor angular velocity Omega [rad/s], Tip angle theta_tw [rad],                |
% |              Cl_alpha [~], Cd_mean [~].                                                     |
% |Output      : Tc, Hc, Yc, Qc, Pc, the angle of attack [deg] and the rotor inflow ratio       |
% |Note        :                                                                                |
% ===============================================================================================
close all; clear all; clc;

% the following script is only aimed to test the Ndim_Coeff_Articulated_Rotor function

% Input values
Np       = 30;                    % Number of points
V_max    = 80;                    % Maximum speed [m/s]
V_inf    = linspace(0,V_max,Np); 
altitude = 3000;                  % [m]
Lock     = 8;                     % [~]  
f        = 0.7785;                % [m^2]
X        = 0;                     % [rad]

N        = 4;                     % Blade number
W       = 8600;                   % Weight [kg]
R        = 7.79;                  % radius [m]
c        = 0.6;                   % chord [m]
Cl_alpha = 2*pi;
Omega    = 27.856;                % [rad/s]
theta_tw = -8;                    % [deg]
Cd_mean  = 0.0121;

sigma    = N*c/(pi*R);            % Rotor solidity
theta_tw = convang(theta_tw,'deg','rad');

for i = 1:length(V_inf)
[Tc(i),Hc(i),Yc(i),Qc(i),Pc(i),alfa(i),lambda(i)] = Ndim_Coeff_Articulated_Rotor(V_inf(i),altitude,Lock,f,X,sigma,W,R,Omega,theta_tw,Cl_alpha,Cd_mean);
end
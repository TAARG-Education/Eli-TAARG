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
% |Name        : main.m                                                                         |
% |Author      : Marco Artiano & Luca Angelino                                                  |
% |              University of Naples Federico II.                                              |
% |Version     : 1.0.1                                                                          |
% |Date        : 24/01/21                                                                       |
% |Modified    : 24/01/21                                                                       |
% |Description : Ndim_Coeff_Articulated_Rotor.m function testing                                |
% |Reference   : 'Lezioni di AERODINAMICA DELL'ALA ROTANTE a.a. 2019-2020 - Renato Tognaccini'  |
% |Input       : Velocity with respect to the air (V_inf), Altitude(altitude)                   |
% |Output      : Tc, Hc, Yc, Qc, Pc and the angle of attack                                     |
% |Note        :                                                                                |
% ===============================================================================================
close all; clear all; clc;

% the following script is only aimed to test the Ndim_Coeff_Articulated_Rotor function

% Input values
V_inf    = 30;      % [m/s]
altitude = 3000;    % [m]
Lock     = 8;       % [~]  
f        = 0.7785;  % [m^2]
X        = 0;       % [rad]

[Tc,Hc,Yc,Qc,Pc,alfa] = Ndim_Coeff_Articulated_Rotor(V_inf,altitude,Lock,f,X);


% % \RVortexInt.m
%  \The function evaluates the dimensionless radial velocity component
%   of the ring vortex by means of integration of Biot-Savart law
%  \Olino Massimiliano, Marino Giuseppe
%  \1.0
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
% ==============================================================================================
% |Name        : RVortexInt.m                                                                  |
% |Author      : Olino Massimiliano, Marino Giuseppe                                           |
% |              University of Naples Federico II.                                             |
% |Version     : 1.0                                                                           |
% |Date        : 01/12/2020                                                                    |
% |Modified    : 26/01/2021                                                                    |
% |Description : The function evaluates the dimensionless radial velocity component            |
% |              of the ring vortex by means of integration of Biot-Savart law                 |
% |Reference   : McCormick, B.W.,(1967), Aerodynamics of V/STOL Flight, Academic Press.        |
% |              Tognaccini, R., (2020), Lezioni di AERODINAMICA DELL'ALA ROTANTE.             |
% |Input       : ch (chord), D14, D34 (The diameters at 1/4 and 3/4 of the shroud)             |
% |Output      : fx (the dimensionless radial velocity component of the ring vortex )          |
% |Note        :                                                                               |
% ==============================================================================================
clear all; close all; clc;

ch=2.40;                        % Chord of the duct
D14=1.50;                       % Diameter at 1/4 of the duct
D34=1.00;                       % Diameter at 3/4 of the duct

f=Feature_RVortexInt(ch,D14,D34);

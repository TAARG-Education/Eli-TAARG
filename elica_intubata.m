%% \elica_intubata.m
%  \brief:  A function that calculates total thrust of a ducted propeller. 
%   It generates a vector with total thrust and shroud thrust as output.
%  \author: Claudio Mirabella, Christian Salzano
%  \version: 1.04
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
% |Name        : elica_intubata.m                                                                  |
% |Author      : Claudio Mirabella, Christian Salzano                                    |
% |              University of Naples Federico II.                                       |
% |Version     : 1.04                                                                    |
% |Date        : 25/11/2020                                                              |
% |Modified    : 11/01/2021                                                              |
% |Description : A function that calculates total thrust of a ducted propeller.          | 
% |              It generates a vector with total thrust and shroud thrust as output.  |     
% |Reference   : Lezioni di Aerodinamica dell'Ala Rotante                                |                                                         
% |Input       :                                                                         |
% |              (rho) = Density                                                         |
% |              (Gamma) = Ring vortex circulation associated with the shroud            |  
% |              (Vinf) = Stream velocity                                                |
% |              (Dquarter) = Shroud diameter at c/4                                     |
% |              (R) = Shroud radius                                                     |
% |              (TR) = Isolated rotor thrust                                            |
% |Output      : (T) = Total thrust generated                                            |
% |              (TS) = Thrust generated by the shroud                                   |
% |Note        :                                                                         |
% ==============================================================================================

function [T, TS] = elica_intubata(rho, Gamma, Vinf, Dquarter, R, TR)
% ---------------------------------------------------------------------
% Propeller's induced velocity calculations
%  eq 4.24
w0 = .5*(-Vinf + sqrt(Vinf^2 + 2*TR/(rho*pi*R^2)));
% ---------------------------------------------------------------------
% Streamtube's radius rt calculated via the continuity
%  eq 4.20
A = pi*R^2;
const=0;
const = (Vinf + w0)*A;
rt = sqrt(const/(Vinf*pi));
% ---------------------------------------------------------------------
% Radial velocity component induced by the rotor at c/4
%  eq 4.23
cquarter=1/4;
viRquarter = -.5*rt*w0*R^2/((R^2+cquarter^2)^1.5);
% ---------------------------------------------------------------------
% Thrust component due to the shroud
TS = -rho*viRquarter*Gamma*pi*Dquarter;
% ---------------------------------------------------------------------
% Total thrust
T = TR + TS;
% ---------------------------------------------------------------------
end

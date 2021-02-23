%% \main.m
%  \brief:  It is a test to validate ducted_prop_thrust.m
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
% |Name        : main.m                                                                  |
% |Author      : Claudio Mirabella, Christian Salzano                                    |
% |              University of Naples Federico II.                                       |
% |Version     : 1.04                                                                    |
% |Date        : 25/11/2020                                                              |
% |Modified    : 18/02/2021                                                              |
% |Description : It is a test to validate elica_intubata.m                               |     
% |Reference   : Lezioni di Aerodinamica dell'Ala Rotante                                |                                                         
% |Input       :                                                                         |
% |              (rho) = 1.225 [kg/m^3]                                                  |
% |              (ch) = 2.4 [m]                                                          | 
% |              (Vinf) = 35 [m/s]                                                       |
% |              (D14) = 1.5 [m]                                                         |
% |              (D34) = 1.5 [m]                                                         |
% |              (theta) = 0 [deg]                                                       |
% |              (R) = .75 [m]                                                           |
% |              (TR) = 500 [N]                                                          |
% |Output      : (T) = 566.3129 [N]                                                      |
% |              (TS) = 66.3129 [N]                                                      |
% |Note        :                                                                         |
% ==============================================================================================

rho = 1.225;
ch = 2.40;
Vinf = 35;
D14 = 1.50;
D34 = 1.50;
theta = 0;                   
R = .75;                         
TR = 500;                       

[T,TS] = ducted_prop_thrust(rho, ch, Vinf, D14, D34, theta, R, TR)

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
% |Modified    : 11/01/2021                                                              |
% |Description : It is a test to validate elica_intubata.m                               |     
% |Reference   : Lezioni di Aerodinamica dell'Ala Rotante                                |                                                         
% |Input       :                                                                         |
% |              (rho) = 1.225 [kg/m^3]                                                  |
% |              (Gamma) = 30 [m^2/s]                                                    |  
% |              (Vinf) = 35 [m/s]                                                       |
% |              (Dquarter) = 2 [m]                                                      |
% |              (R) = .9 [m]                                                            |
% |              (TR) = 500 [N]                                                          |
% |Output      : (T) = 729.6644 [N]                                                      |
% |              (TS) = 229.6644 [N]                                                     |
% |Note        :                                                                         |
% ==============================================================================================

Dquarter = 2;                   
Gamma = 30;
Vinf = 35;                      
rho = 1.225;                    
R = .9;                         
TR = 500;                       

[T,TS] = ducted_prop_thrust(rho,Gamma,Vinf,Dquarter,R,TR)

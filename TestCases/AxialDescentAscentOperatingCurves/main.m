%% \Axial_Descent_Ascent_Operating_Curves_Rotor.m
%  \brief: the function plots w(V_infty) and P(V_infty) curves according to
%  Impulsive theory.
%  aerodynamic model
%  \author: Colledà Moreno, Veneruso Salvatore
%  \version: 1.00
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
% |Name        : Axial_Descent_Ascent_Operating_Curves_Rotor.m                           |
% |Author      : Colledà Moreno, Veneruso Salvatore                                      |
% |              University of Naples Federico II.                                       |
% |Version     : 1.00                                                                    |
% |Date        : 23/12/2020                                                              |
% |Modified    : 23/12/2020                                                              |
% |Description : the function plots w(V_infty) and P(V_infty) curves according to
%                Rotor Simply Impulsive theory.                                          |     
% |Reference   : Renato Tognaccini. Appunti Aerodinamica dell'ala rotante.
%                Università degli studi di Napoli Federico II. a.a.2020/2021             |
% |Input       : * the inputs must be Mass of rotorcraft and radius of rotor             | 
% |Output      : w(V_infty) and P(V_infty) plots                                         |
% |Note        :                                                                         |
% ==============================================================================================
clc; clear; close all;
M = 5000; %kg;
R = 7;    %m;
[P,I] = Axial_Descent_Ascent_Operating_Curves_Rotor(M,R);



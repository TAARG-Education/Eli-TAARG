%% \Opti_Prop_TEST.m
%  \brief: Test function
%  \authors : Francesco Gervasio, Nicola Russo
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
% |Name        : Opti_Prop.m                                                   |
% |Authors     : Francesco Gervasio, Nicola Russo                              |
% |              University of Naples Federico II.                             |
% |Version     : 1.04                                                          |
% |Date        : 28/11/2020                                                    |
% |Modified    : 23/12/2020                                                    |
% |Description : this function provides the axial and the rotational           |
% inductions of the optimal propeller as well as the thrust and power          |
% coefficient distributions vs the adimensional radius                         |
% for a given thrust or power cofficients.                                     |     
% |Reference   :                                                               |
% |                                                                            |                                                         |
% |Input     : GEOMETRY
% |            N_blade = Number of blades         [ ]                          |  
% |            R_hub   = Hub radius percentage    [%]                          |
% |            n_rpm   = Revolutions per minutes  [rpm]                        |
% |              ATMO                                                          |
% |              h       = Height               [m]                            |
% |              AERO                                                          |
% |              V_inf   = Asymptotic speed     [m/s]                          |
% |              Ct       = Thrust coefficients                                |
% |              Cp       = Power coefficients                                 |

% |Note        :                                                               |    
% ==============================================================================================
clear all; close all; clc;

%% INPUT
N_blade = 2;        %[ ]
R_hub   = 0.2;      %[%]
R_tip   = 1.2;      %[m]
n_rpm   = 2000;     %[rpm]
V_inf   = 45 ;      %[m/s]
Ct = 0.03;          %[]
Cp =0.0183;         %[]
h = 4000;           %[m]

[r_adim_T,a_chi_T,a_first_chi_T,dCt,w0_T,k_T,Ct_new,dCp] = Opti_prop_T(N_blade,R_hub,R_tip,n_rpm,V_inf,Ct,h); 
[r_adim_P,a_chi_P,a_first_chi_P,dCp,w0_P,k_P,Cp_new,dCt] = Opti_prop_P(N_blade,R_hub,R_tip,n_rpm,V_inf,Cp,h); 



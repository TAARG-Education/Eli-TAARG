%% \Opti_Prop_TEST.m
%  \brief: Test function
%  \authors : Francesco Gervasio, Nicola Russo
%  \version: 1.03
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
% |Version     : 1.03                                                         |
% |Date        : 28/11/2020                                                    |
% |Modified    : 20/12/2020                                                     |
% |Description : this function designs the geometry of an optimum propeller
%  at given thrust or power.                                                   |     
% |Reference   :                                                               |
% |                                                                            |                                                         |
% |Input       : GEOMETRY
% |              N_blade = Number of blades [ ]                                |  
% |              R_hub   = Hub radius percentage  [%]                          |
% |              n_rpm   = Revolutions per minutes  [rpm]                      |
% |              ATMO                                                          |
% |              h       = Height               [m]                            |
% |              AERO                                                          |
% |              V_inf   = Asymptotic speed     [m/s]                          |
% |              T       = Thrust               [N]                            |
% |              P       = Power                [W]                            |
% |              Aerodynamic characteristics at a given section                |
% |                                                                            |
% |                                                                            |
% |Output      : (Theta) = pitch angle distribution                            |
% |                                                                            |
% |                                                                            |
% |Note        :                                                               |
% ==============================================================================================
clear all; close all; clc;

%% INPUT
N_blade = 2;        %[ ]
R_hub   = 0.2;     %[%]
n_rpm   = 2000;     %[rpm]
V_inf   = 45 ;      %[m/s]
Ct = 0.05;           %[]
Cp =0.060;           %[]
h = 4000;           %[m]

[r_adim_T_new,a_chi,a_first_chi,dCt_N,W_0,k,Ct_corr] = Opti_prop_T_Newton(N_blade,R_hub,n_rpm,V_inf,Ct,h);
[r_adim_T,a_chi,a_first_chi,dCt,w0] = Opti_prop_T_corr(N_blade,R_hub,n_rpm,V_inf,Ct,h); 
error_perc = abs((Ct_corr - Ct))/Ct*100

% figure
% plot(r_adim_T_new,dCt_N,'-k');
% grid on;
% hold on ;
% plot(r_adim_T,dCt,'ok');
% xlabel('$\bar{r}$','interpreter','latex');
% ylabel('$\frac{dC_{T}}{d\bar{r}}$','interpreter','latex');
% title('Thrust coefficient distribution $\frac{dC_{T}}{d\bar{r}}$ vs $\bar{r}$','interpreter','latex');
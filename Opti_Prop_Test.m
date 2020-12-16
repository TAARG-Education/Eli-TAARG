%% \Opti_Prop_TEST.m
%  \brief: Test function
%  \authors : Francesco Gervasio, Nicola Russo
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
% ==============================================================================================
% |Name        : Opti_Prop.m                                                   |
% |Authors     : Francesco Gervasio, Nicola Russo                              |
% |              University of Naples Federico II.                             |
% |Version     : 1.01                                                          |
% |Date        : 28/11/2020                                                    |
% |Modified    : 4/12/2020                                                     |
% |Description : this function designs the geometry of an optimum propeller
%  at given thrust or power.                                                   |     
% |Reference   :                                                               |
% |                                                                            |                                                         |
% |Input       : GEOMETRY
% |              R_prop  = Propeller radius [m]                                |
% |              N_blade = Number of blades [ ]                                |  
% |              R_hub   = Hub radius       [m]                                |
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
R_tip  = 0.95;      %[m] 
N_blade = 2;        %[ ]
R_hub   = 0.14;     %[m]
n_rpm   = 2500;     %[rpm]
V_inf   = 45 ;      %[m/s]
T = 964;            %[N]
P = 53579.12;       %[kW]
h = 4000;           %[m]

[r_adim_T,theta_T,chord_T,dT_corr,w_0_T] = Opti_prop_T(R_tip,N_blade,R_hub,n_rpm,V_inf,T,h);
[r_adim_P,theta_P,chord_P,dP_corr,w_0_P] = Opti_prop_P(R_tip,N_blade,R_hub,n_rpm,V_inf,P,h);

figure
plot(r_adim_T,theta_T,'-k');
hold on;
plot(r_adim_P,theta_P,'ok','markersize',3);
xlabel('$\bar{r}$','interpreter','latex');
ylabel('$\theta$ (rad)','interpreter','latex');
title('Theta angle $ \theta$  vs $\bar{r}$','interpreter','latex');
grid on;

figure
plot(r_adim_T,chord_T,'-k');
hold on;
plot(r_adim_P,chord_P,'ok','markersize',3);
xlabel('$\bar{r}$','interpreter','latex');
ylabel('$c(\bar{r})$','interpreter','latex');
title('Chord distribution c vs $\bar{r}$','interpreter','latex');
%% \Opti_Prop_TEST.m
%  \brief: Test function
%  \authors : Francesco Gervasio, Nicola Russo
%  \version: 1.08
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
% |Version     : 1.08                                                          |
% |Date        : 28/11/2020                                                    |
% |Modified    : 29/01/2021                                                    |
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
R_hub   = 0.15;     %[%]
R_tip   = 0.9;      %[m]
n_rpm   = 2500;     %[rpm]
V_inf   = 58.33 ;   %[m/s]
Ct = 0.0740;        %[]
h = 4510;           %[m]

r = linspace(R_hub*R_tip,R_tip,1000)/R_tip;
[r_adim_T,chi_T,a_corr_T,a_first_corr_T,dCt_dradim_T,dCp_dradim_T,Cp] = Opti_prop_T(N_blade,R_hub,R_tip,n_rpm,V_inf,Ct,h); 
[r_adim_P,chi_P,a_corr_P,a_first_corr_P,dCp_dradim_P,dCt_dradim_P,Ct] = Opti_prop_P(N_blade,R_hub,R_tip,n_rpm,V_inf,Cp,h); 


%% Plot
figure
plot(chi_T,a_corr_T,'-k')
hold on
plot(chi_T,a_first_corr_T,'--k')
hold on;
plot(chi_P,a_corr_P,'ob','markersize',2)
hold on
plot(chi_P,a_first_corr_P,'or','markersize',2)
xlabel('$\chi$','interpreter','latex');
legend('$a_{T}$','$a^{''}_{T}$','$a_{P}$','$a^{''}_{P}$','interpreter','latex')
;grid on;
title('Axial induction $ a$ and rotational induction $a^{''}$ vs $\chi$','interpreter','latex');

figure
plot(r_adim_T,dCt_dradim_T,'-k');
hold on;
plot(r_adim_P,dCt_dradim_P,'ob','markersize',2);
grid on;
hold on;
legend('Opti\_Prop\_T','Opti\_Prop\_P','interpreter','latex')
xlabel('$\bar{r}$','interpreter','latex');
ylabel('$\frac{dC_{T}}{d\bar{r}}$','interpreter','latex');
title('Thrust coefficient distribution $\frac{dC_{T}}{d\bar{r}}$ vs $\bar{r}$','interpreter','latex');


figure
plot(r_adim_T,dCp_dradim_T,'-k');
hold on;
plot(r_adim_P,dCp_dradim_P,'ob','markersize',2);
grid on;
hold on;
xlabel('$\bar{r}$','interpreter','latex');
ylabel('$\frac{dC_{P}}{d\bar{r}}$','interpreter','latex');
legend('Opti\_Prop\_T','Opti\_Prop\_P','interpreter','latex')
title('Power coefficient distribution $\frac{dC_{P}}{d\bar{r}}$ vs $\bar{r}$','interpreter','latex');
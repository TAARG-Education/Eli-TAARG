%% \Opti_Prop.m
%  \brief: this function designs the geometry of an optimum propeller
%  at given thrust or power.
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
% |              P       = Power                [kW]                           |
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
R_prop  = 0.454;    %[m] 
N_blade = 2;        %[ ]
R_hub   = 0.07;     %[m]
A = pi * R_prop^2;  %[m]^2
D_prop = R_prop*2;  %[m]
n_rpm   = 2000;     %[rpm]
V_inf   = 27.7 ;    %[m/s]
T = 150000;         %[N]
P = 500;            %[kW]
h = 2000;           %[m]

%% CONVERSIONS
%{
From rpm to rps and rad/s
%}

n_rps = n_rpm/60; 
n_rads = convangvel(n_rpm,'rpm','rad/s');

%% ATMO
%{ 
The aerospace toolbox is required in order to calculate the ISA parameters.
It can be downloaded from the Mathworks website.
%}
[T a_chi p rho] = atmosisa(h);

%% Ct_Cp
%{
Given the required thrust or the required power, the related coefficients
are computed.
%}
Ct = T/(rho*n_rps^2*D_prop^4);
%Cp =P/(rho*n_rps^3*D_prop^5);

%% Step 1 
%{
The first attempt value of the axial induction w_0 has been set to w_0 =
a*V_inf where a is the axial induction a = w_0/V_inf. 
%}
syms w0;
eqn1 = w0 > 0;
eqn2 = T == 2*rho*A*(V_inf+w0)*w0;
eqn = [eqn1 eqn2];
S = solve(eqn,w0,'ReturnConditions',true);
w_0 = double(S.w0);

%% Step 2
%{
The inflow angle distribution in already determined. A plot of this
parameter has been provided as well.
%}
r = linspace(R_hub,R_prop,1000);
r_adim = r/R_prop;
phi = atan((V_inf+w_0)./(n_rads*r))';
phi_adim = atan((V_inf+w_0)./(n_rads*r_adim))';
figure
plot(r_adim,phi_adim,'-k');
xlabel('$\bar{r}$','interpreter','latex');
ylabel('$\phi$ (rad)','interpreter','latex');
title('Inflow angle $ \phi$  vs $\bar{r}$','interpreter','latex');
grid on;

%% Step 3
%{
In this step a calculation of:
- chi = (OMEGA*r)/V_inf (an alternative non-dimensional radius)
- a   = axial induction as function of both chi and adim radius r_adim
- a'  = rotational induction as function of both chi and adim radius r_adim
has been provided and plotted.
%}
chi = (n_rads.*r)/V_inf;
a_chi = w_0/V_inf*((chi.^2)./((1+w_0/V_inf)^2+(chi.^2)));
a_first_chi = w_0/V_inf*((1+w_0/V_inf)./((1+w_0/V_inf)^2+chi.^2));
a_radim = w_0/V_inf*(((n_rads*r_adim*R_prop/V_inf).^2)./((1+w_0/V_inf)^2+((n_rads*r_adim*R_prop/V_inf).^2)));
a_first_radim = w_0/V_inf*((1+w_0/V_inf)./((1+w_0/V_inf)^2+(n_rads*r_adim*R_prop/V_inf).^2));
figure
plot(chi,a_chi,'-.k')
hold on
plot(chi,a_first_chi,'-k')
xlabel('$\chi$','interpreter','latex');
legend('$a$','$a^{''}$','interpreter','latex');
title('Axial induction $ a$ and rotational induction $a^{''}$ vs $\bar{r}$','interpreter','latex');

%{
rotational and flight speed velocity have been composed at any given
section in order to obtain the true velocity. A plot has been provided as...
well 
%}
V_eff = (V_inf+w_0)./sin(phi);
figure
plot(r_adim,V_eff,'-k');
xlabel('$\bar{r}$','interpreter','latex');
ylabel('$V_{eff}$ (m/s)','interpreter','latex');
title('Velocity $V_{eff}$  vs $\bar{r}$','interpreter','latex');
grid on;
hold on;


%{
The Prandtl correction function for finite blade number has been
implemented. Lambda is  = chi^-1. According to the momentum theory: ...
w_j = 2*w_0;
%}
lambda = (n_rads*R_prop/V_inf)^-1;
w_j = 2*w_0;
F   = (2/pi)*acos(exp((N_blade/(2*lambda))*(r_adim-1)));
figure
plot(r_adim,F,'-k')
xlabel('$\bar{r}$','interpreter','latex');
ylabel('$F$','interpreter','latex');
title('Prandtl correction function for finite blade number vs $\bar{r}$',...
    'interpreter','latex');
grid on;

%{
The non-dimensional aerodynamic load gamma_adim_opti, already scaled with the Prandtl
correction function, can be plotted as well. The non dimensional Gamma function can be
visualized in the aforementioned plot as ratio.
%} 
GAMMA = (4*pi.*F.*r_adim.^2.*a_first_radim*n_rads*R_prop^2)/N_blade;
gamma_adim_opti = (N_blade.*GAMMA*n_rads)/(4*pi*V_inf*w_0);
figure
plot(chi,gamma_adim_opti,'-k')
grid on;
xlabel('$\chi$','interpreter','latex');
legend('$\frac{N \Gamma \Omega}{4\pi w_{0}*v_{inf}}$','interpreter','latex');


%{
At each step, dT/dr_adim or dP/dr_adim are calculated. The optimal
distribution is supposed to be linear. These values are lately integrated
along the non dimensional radius in order to calculate a step T or P value that
must be confronted with the design T or P value.
%}


sigma_cl = (1/pi)*((N_blade.*GAMMA)./(V_eff'.*r_adim.*R_prop));

% %% Step 3
% cl =  [0.8 0.8 0.8 0.6 0.2 0];
% cd =  [100 100 100 100 100 0]./10000;
% r  =  [0.2 0.4 0.6 0.8 0.9 1];
% mu = 1e-5;
% c = 1;
% Re =  (rho.*V_eff*c)/mu;
% RE = [Re(200) Re(400) Re(600) Re(800) Re(900) Re(1000)];
% alpha = [2.8 2.79 2.77 1 -2.45 -4.19]/57.63;
% alpha_int = @(t)interp1(r,alpha,t);
% cl_int = @(t)interp1(r,cl,t);
% cd_int = @(t)interp1(r,cd,t);
% 
% t = linspace(0.20,1,1000);
% alpha_vec = alpha_int(t);
% cl_vec = cl_int(t);
% cd_vec = cd_int(t);
% theta = phi(201:1000)+alpha_vec(201:1000);
% 
% %% Calcolo perdite viscose
% 
% lambda_first = cl_vec.*cos(phi')-cd_vec.*sin(phi');
% 
% dT_dr_adim = (0.5*sigma_cl.*r_adim.^2).*lambda_first;
% T_star = trapz(dT_dr_adim(2:1000));
% figure
% plot(r_adim,dT_dr_adim);



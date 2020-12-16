function [r_adim,theta,chord,dP_corr,w_0] = Opti_prop_P(R_tip,N_blade,R_hub,n_rpm,V_inf,P,h)
%% \Opti_Prop.m
%  \brief: this function designs the geometry of an optimum propeller
%  at given Power.
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

%% Geometry
A = pi * R_tip^2;   %[m]^2


%% CONVERSIONS
%{
From rpm to rps and rad/s
%}
n_rads = convangvel(n_rpm,'rpm','rad/s');

%% ATMO
%{ 
The aerospace toolbox is required in order to calculate the ISA parameters.
It can be downloaded from the Mathworks website.
%}
[Temp a p rho] = atmosisa(h);

%% Step 1 
%{
The first attempt value of the axial induction w_0 has been set to w_0 =
a*V_inf where a is the axial induction a = w_0/V_inf. 
%}
syms w0;
eqn1 = w0 > 0;
eqn2 = P == 2*rho*A*((V_inf+w0).^2)*w0;
eqn = [eqn1 eqn2];
S = solve(eqn,w0,'ReturnConditions',true);
w_0 = double(S.w0);
P_corr = 0;

while (P-P_corr)/P >0.0001
%% Step 2
%{
The inflow angle distribution in already determined. A plot of this
parameter has been provided as well.
%}
r = linspace(R_hub,R_tip,1000);
r_adim = r/R_tip;
phi = atan((V_inf+w_0)./(n_rads*r));
phi_adim = atan((V_inf+w_0)./(n_rads*r_adim))';

%% Step 3.1
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
a_radim = w_0/V_inf*(((n_rads*r_adim*R_tip/V_inf).^2)./((1+w_0/V_inf)^2+((n_rads*r_adim*R_tip/V_inf).^2)));
a_first_radim = w_0/V_inf*((1+w_0/V_inf)./((1+w_0/V_inf)^2+(n_rads*r_adim*R_tip/V_inf).^2));

%{
rotational and flight speed velocity have been composed at any given
section in order to obtain the true velocity. A plot has been provided as...
well 
%}
V_eff = (V_inf+w_0)./sin(phi);

%{
The Prandtl correction function for finite blade number has been
implemented. Lambda is  = chi^-1. According to the momentum theory: ...
w_j = 2*w_0;
%}
lambda = (n_rads*R_tip/V_inf)^-1;
w_j = 2*w_0;
F   = (2/pi)*acos(exp((N_blade/(2*lambda))*(r_adim-1)));

%{
The non-dimensional aerodynamic load gamma_adim_opti, already scaled with...
the Prandtl correction function, can be plotted as well. The non ...
dimensional Gamma function can be visualized in the aforementioned plot...
 as ratio.
%} 

GAMMA = (4*pi.*F.*r_adim.^2.*a_first_radim*n_rads*R_tip^2)/N_blade';
gamma_adim_opti = (N_blade.*GAMMA*n_rads)/(4*pi*V_inf*w_0);

%% Step 3.2
%{
At each step, dT/dr_adim or dP/dr_adim are calculated. The optimal
distribution is supposed to be linear. These values are lately integrated
along the non dimensional radius in order to calculate a step T or P ...
value that must be confronted with the design T or P value.
%}
dr_adim = gradient(r_adim);
dP = N_blade*rho*V_inf*(R_tip^2)*n_rads*r_adim.*(1+a_radim).*GAMMA.*dr_adim;
Ph = trapz(dP);

%% Step 4 (Xfoil Interface function in needed)
%{
In this step the aerodynamics characteristics must be obtained by using the
related xfoil interface function. It is necessary to choose one or more
airfoil along the blade span and perform the standard xfoil procedure to
obtain the cl_{id} and the corresponding cd values. The angle alpha value
must be provided as well as the Reynolds number at the related section
%}
sigma_cl = (1/pi)*((N_blade.*GAMMA)./(V_eff.*r_adim*R_tip));
cl =  [0.8 0.8 0.8 0.6 0.4 ];
cd =  [100 100 100 100 100 ]./10000;
r_1=  [0.2 0.4 0.6 0.8 0.9 ];
mu = 1e-5;
c = 1;
Re =  (rho.*V_eff*c)/mu;
RE = [Re(200) Re(400) Re(600) Re(800) Re(900) Re(1000)];
alpha = [2.8 2.79 2.77 1 -2.45]/57.63;
alpha_int = @(t)interp1(r_1,alpha,t,'spline','extrap');
cl_int = @(t)interp1(r_1,cl,t,'spline','extrap');
cd_int = @(t)interp1(r_1,cd,t,'spline','extrap');
t = linspace(0,1,1000);
alpha_vec = alpha_int(t);
cl_vec = cl_int(t);
cd_vec = cd_int(t);
theta = phi+alpha_vec;

% Viscous losses calculation
lambda_second = cl_vec.*sin(phi)+cd_vec.*cos(phi);

% Propeller solidity calculation
sigma_radim = sigma_cl./cl_vec;
chord = R_tip.*((2*pi.*r_adim.*sigma_radim)/N_blade);

%% Step 5
%{
The thrust/power has been calculating by considering the aforementioned ...
viscous loss; If the while loop condition is not verified, the the cycle is
repeted with an incremented value of the induction velocity.
%}
dP_corr = rho*pi*sigma_radim.*lambda_second.*n_rads.*(r_adim.^2).*(V_eff.^2).*dr_adim;
P_corr = trapz(dP_corr);
w_0 = w_0 +0.0001;
end

% %% Plot
% figure
% plot(r_adim,phi_adim,'-k');
% xlabel('$\bar{r}$','interpreter','latex');
% ylabel('$\phi$ (rad)','interpreter','latex');
% title('Inflow angle $ \phi$  vs $\bar{r}$','interpreter','latex');
% grid on;
% 
% figure
% plot(chi,a_chi,'-.k')
% hold on
% plot(chi,a_first_chi,'-k')
% xlabel('$\chi$','interpreter','latex');
% legend('$a$','$a^{''}$','interpreter','latex');
% title('Axial induction $ a$ and rotational induction $a^{''}$ vs $\bar{r}$','interpreter','latex');
% 
% figure
% plot(r_adim,F,'-k')
% xlabel('$\bar{r}$','interpreter','latex');
% ylabel('$F$','interpreter','latex');
% title('Prandtl correction function for finite blade number vs $\bar{r}$',...
%     'interpreter','latex');
% grid on;
% 
% figure
% plot(chi,gamma_adim_opti,'-k')
% grid on;
% xlabel('$\chi$','interpreter','latex');
% legend('$\frac{N \Gamma \Omega}{4\pi w_{0}*v_{inf}}$','interpreter','latex');
% 
% figure
% plot(r_adim,theta,'-k');
% xlabel('$\bar{r}$','interpreter','latex');
% ylabel('$\theta$ (rad)','interpreter','latex');
% title('Theta angle $ \theta$  vs $\bar{r}$','interpreter','latex');
% grid on;
% 
% figure
% plot(r_adim,chord,'-k');
% grid on;
% xlabel('$\bar{r}$','interpreter','latex');
% ylabel('$c(\bar{r})$','interpreter','latex');
% title('Chord distribution c vs $\bar{r}$','interpreter','latex');
% 
% figure
% plot(r_adim,dP_corr,'-k');
% grid on;
% xlabel('$\bar{r}$','interpreter','latex');
% ylabel('$\frac{dT}{d\bar{r}}$','interpreter','latex');
% title('Power distribution $P$ vs $\bar{r}$','interpreter','latex');

DATA(:,1) = r_adim';
DATA(:,2) = theta';
DATA(:,3) = chord';
DATA(:,4) = dP_corr';
filename = ['Data_Opti_Prop_P.txt'];
fid = fopen(filename, 'wt');
fprintf(fid, '%s\t%s\t%s\t%s\n', '   r_adim','theta(rad)','chord(r_adim)','dP_corr(r_adim)');  % header
fclose(fid);
dlmwrite(filename,DATA,'delimiter','\t','precision',['%10.',num2str(4),'f'],'-append');
end

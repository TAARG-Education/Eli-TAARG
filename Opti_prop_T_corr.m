function [r_adim,a_chi,a_first_chi,dCt_corr,w_0] = Opti_prop_T_corr(N_blade,R_hub,n_rpm,V_inf,Ct,h)
%% \Opti_Prop.m
%  \brief: this function designs the geometry of an optimum propeller
%  at given thrust.
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
% |              Ct       = Thrust coefficient  [N]                            |
% |              Cp       = Power  coefficient  [W]                            |
% |              Aerodynamic characteristics at a given section                |
% |                                                                            |
% |                                                                            |
% |Output      : (Theta) = pitch angle distribution                            |
% |                                                                            |
% |                                                                            |
% |Note        :                                                               |
% ==============================================================================================

%% INPUT
%{
R_tip  =    [m] 
N_blade =   [ ]
R_hub   =   [%]
n_rpm   =   [rpm]
V_inf   =   [m/s]
T =         [N]
P =         [kW]
h =         [m]
%}

%% CONVERSIONS
%{
From rpm to rps and rad/s
%}
n_rads = convangvel(n_rpm,'rpm','rad/s');
n_rps = n_rpm/60;

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
T = Ct*rho*(n_rps^2)*(16);
syms w0;
eqn1 = w0 > 0;
eqn2 = T == 2*rho*pi*(V_inf+w0)*w0;
eqn = [eqn1 eqn2];
S = solve(eqn,w0,'ReturnConditions',true);
w_0 = double(S.w0);
Ct_corr = 0;
r_adim = linspace(R_hub,1,1000);
chi = (n_rads.*r_adim)/V_inf;


while (Ct-Ct_corr)/Ct >0.0001
%% Step 3.1
%{
In this step a calculation of:
- chi = (OMEGA*r)/V_inf (an alternative non-dimensional radius)
- a   = axial induction as function of both chi and adim radius r_adim
- a'  = rotational induction as function of both chi and adim radius r_adim
has been provided and plotted.
%}

a_chi = w_0/V_inf*((chi.^2)./((1+w_0/V_inf)^2+(chi.^2)));
a_first_chi = w_0/V_inf*((1+w_0/V_inf)./((1+w_0/V_inf)^2+chi.^2));
a_radim = w_0/V_inf*(((n_rads*r_adim/V_inf).^2)./((1+w_0/V_inf)^2+((n_rads*r_adim/V_inf).^2)));
a_first_radim = w_0/V_inf*((1+w_0/V_inf)./((1+w_0/V_inf)^2+(n_rads*r_adim/V_inf).^2));

%{
The Prandtl correction function for finite blade number has been
implemented. Lambda is  = chi^-1. According to the momentum theory: ...
w_j = 2*w_0;
%}
lambda = (n_rads/V_inf)^-1;
w_j = 2*w_0;
F   = (2/pi)*acos(exp((N_blade/(2*lambda))*(r_adim-1)));

%{
The non-dimensional aerodynamic load gamma_adim_opti, already scaled with...
the Prandtl correction function, can be plotted as well. The non ...
dimensional Gamma function can be visualized in the aforementioned plot...
 as ratio.
%} 

GAMMA = (4*pi.*F.*r_adim.^2.*a_first_radim*n_rads)/N_blade';
gamma_adim_opti = (N_blade.*GAMMA*n_rads)/(4*pi*V_inf*w_0);


%% Step 3.2

%{
At each step, dT/dr_adim or dP/dr_adim are calculated. The optimal
distribution is supposed to be linear. These values are lately integrated
along the non dimensional radius in order to calculate a step T or P ...
value that must be confronted with the design T or P value.
%}
dr_adim = gradient(r_adim);
dT = N_blade*rho*n_rads*r_adim.*(1-a_first_radim).*GAMMA.*dr_adim;
dCt_corr = dT/(rho*(n_rps^2)*(2^4));
Ct_corr = trapz(dCt_corr);
w_0 = w_0 +0.001;


end
%% Plot
figure
plot(chi,a_chi,'-.k')
hold on
plot(chi,a_first_chi,'-k')
xlabel('$\chi$','interpreter','latex');
legend('$a$','$a^{''}$','interpreter','latex');
grid on;
title('Axial induction $ a$ and rotational induction $a^{''}$ vs $\bar{r}$','interpreter','latex');

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
figure
plot(r_adim,dCt_corr,'-k');
grid on;
xlabel('$\bar{r}$','interpreter','latex');
ylabel('$\frac{dC_{T}}{d\bar{r}}$','interpreter','latex');
title('Thrust coefficient distribution $\frac{dC_{T}}{d\bar{r}}$ vs $\bar{r}$','interpreter','latex');

DATA(:,1) = r_adim';
DATA(:,2) = chi';
DATA(:,3) = a_chi';
DATA(:,4) = a_first_chi';
DATA(:,5) = dCt_corr';
filename = ['Data_Opti_Prop_T.txt'];
fid = fopen(filename, 'wt');
fprintf(fid, '%s\t%s\t%s\t%s\t%s\n', '  r_adim','  chi','      a(chi)','  a''(chi)','  dCt_corr(r_adim)');  % header
fclose(fid);
dlmwrite(filename,DATA,'delimiter','\t','precision',['%10.',num2str(6),'f'],'-append');
end



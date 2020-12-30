%% \Opti_Prop_TEST.m
%  \brief: Test function
%  \authors : Francesco Gervasio, Nicola Russo
%  \version: 1.05
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
% |Version     : 1.05                                                          |
% |Date        : 28/11/2020                                                    |
% |Modified    : 30/12/2020                                                    |
% |Description : this function provides the axial and the rotational           |
% inductions of the optimal propeller as well as the thrust and power          |
% coefficient distributions vs the adimensional radius                         |
% for a given thrust cofficients.                                              |     
% |Reference   :                                                               |
% |                                                                            |                                                         |
% |Input     : GEOMETRY
% |            N_blade = Number of blades         [ ]                          |  
% |            R_hub   = Hub radius percentage    [%]                          |
% |            R_tip   = Tip radius               [m]                          |
% |            n_rpm   = Revolutions per minutes  [rpm]                        |
% |              ATMO                                                          |
% |              h       = Height               [m]                            |
% |              AERO                                                          |
% |              V_inf   = Asymptotic speed     [m/s]                          |
% |              Ct       = Thrust coefficients                                |
% |Note        :                                                               |
% ==============================================================================================

%% INPUT
%{
N_blade =   [ ]
R_hub   =   [%]
R_tip   =   [m] 
n_rpm   =   [rpm]
V_inf   =   [m/s]
T =         [N]
P =         [kW]
h =         [m]
%}

function [r_adim,a,a_first,dCt_new,w_0,k,Ct_new,dCp] = Opti_prop_T(N_blade,R_hub,R_tip,n_rpm,V_inf,Ct,h)
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

%% Radius and adimensional values

r = linspace(R_hub*R_tip,R_tip,1000);
D = 2*R_tip;
chi = (n_rads.*r)/V_inf;
r_adim = r./R_tip;

%% Step 1 First attempt values
%{
The two first attempt values of the axial induction w_0 and w_1 have been 
set. w_0 = a*V_inf where a is the axial induction a = w_0/V_inf while 
w_1 = 2*w_0. Two different first attempt values are necessary in order to 
initialize the false position method.
%}

T    = Ct*rho*(n_rps^2)*(D^4);
syms w0;
eqn1 = w0 > 0;
eqn2 = T == 2*rho*pi*(V_inf+w0)*w0;
eqn  = [eqn1 eqn2];
S    = solve(eqn,w0,'ReturnConditions',true);
w_0  = double(S.w0);
w_1  = 2*w_0;
T_1  = 2*rho*pi*(V_inf+w_1)*w_1;
Ct_1 = T_1/(rho*(n_rps^2)*(D^4));
Ct_0 = 0;
Ct_new = 0;

%% Errors
tao = 1*10^-8;
error_0 = Ct-Ct_0;
error_1 = Ct-Ct_1;

%% False position method initialization
%{
The false position method has been chosen in ordert to implement the
iterative cycle. The first value of the axial induction speed, has been
calculated according to the classical method formulation by using the first
attempt axial induction speed values.
%}
w = (w_1*error_0-w_0*error_1)/(error_0-error_1);
k = 0; %cycle counter


while abs(Ct-Ct_new) > tao
%% Step 1
%{
In this step a calculation of:
- chi = (OMEGA*r)/V_inf (an alternative non-dimensional radius)
- a   = axial induction as function of chi 
- a'  = rotational induction as function of chi 
has been provided and plotted.
%}
a       = w/V_inf*((chi.^2)./((1+w/V_inf)^2+(chi.^2)));
a_first = w/V_inf*((1+w/V_inf)./((1+w/V_inf)^2+chi.^2));

%{
The Prandtl correction function for finite blade number has been
implemented. Lambda is  = chi^-1. According to the momentum theory: ...
w_j = 2*w;
%}
lambda = (n_rads*R_tip/V_inf)^-1;
F      = (2/pi)*acos(exp((N_blade/(2*lambda))*(r_adim-1)));

%{
The non-dimensional aerodynamic optimal load Gamma, already scaled with
the Prandtl correction function is implemented.
%} 
GAMMA = (4*pi.*F.*(r_adim*R_tip).^2.*a_first*n_rads)/N_blade';

%% Step 2
%{
At each step, dT/dr_adim or dP/dr_adim are calculated. These values are
lately integrated along the non dimensional radius in order to calculate a 
step T or P value that must be confronted with the design T or P value.
%}
dr_adim   = gradient(r_adim);
dT        = N_blade*rho*n_rads*(R_tip^2)*r_adim.*(1-a_first).*GAMMA.*dr_adim;
dCt_new   = dT/(rho*(n_rps^2)*(D^4));
Ct_new    = trapz(dCt_new);
error_new = Ct_new - Ct;

%% Step 3 - FALSE POSITION METHOD
%{
In this step, the induction speed values and errors are updated at each 
cycle and then the new induction speed value is calculated in order to 
reiterate the calculation until the while loop exit condition is verified.
%}
w_0     = w_1;
w_1     = w;
error_0 = error_1;
error_1 = error_new;
w       = (w_1*error_0-w_0*error_1)/(error_0-error_1);
k       = k+1;

end

dCp = dT*(V_inf+w)/(rho*(n_rps^3)*(D^5));
error_perc_Ct = abs((error_new))/Ct*100

%% Plot
figure
plot(chi,a,'-.k')
hold on
plot(chi,a_first,'-k')
xlabel('$\chi$','interpreter','latex');
legend('$a$','$a^{''}$','interpreter','latex');
grid on;
title('Axial induction $ a$ and rotational induction $a^{''}$ vs $\chi$','interpreter','latex');

figure
plot(r_adim,dCt_new,'-k');
grid on;
xlabel('$\bar{r}$','interpreter','latex');
ylabel('$\frac{dC_{T}}{d\bar{r}}$','interpreter','latex');
title('Thrust coefficient distribution $\frac{dC_{T}}{d\bar{r}}$ vs $\bar{r}$','interpreter','latex');

figure
plot(r_adim,dCp,'-k');
grid on;
xlabel('$\bar{r}$','interpreter','latex');
ylabel('$\frac{dC_{P}}{d\bar{r}}$','interpreter','latex');
title('Power coefficient distribution $\frac{dC_{P}}{d\bar{r}}$ vs $\bar{r}$','interpreter','latex');

%% OUTPUT FILE
%{
An output .txt file is created
%}

DATA(:,1) = r_adim';
DATA(:,2) = chi';
DATA(:,3) = a';
DATA(:,4) = a_first';
DATA(:,5) = dCt_new';
DATA(:,6) = dCp';
filename = ['Data_Opti_Prop_T.txt'];
fid = fopen(filename, 'wt');
fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\n', '  r_adim','  chi','      a(chi)','  a''(chi)','  dCt_corr(r_adim)','  dCp_corr(r_adim)');  % header
fclose(fid);
dlmwrite(filename,DATA,'delimiter','\t','precision',['%10.',num2str(6),'f'],'-append');
end
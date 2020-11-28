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

%% CONVERSIONI
n_rps = n_rpm/60;
n_rads = convangvel(n_rpm,'rpm','rad/s');

%% ATMO
[T a_chi p rho] = atmosisa(h);

%% Ct_Cp
Ct = T/(rho*n_rps^2*D_prop^4);
%Cp =P/(rho*n_rps^3*D_prop^5);

%% Step 1
syms w0;
eqn1 = w0 > 0;
eqn2 = T == 2*rho*A*(V_inf+w0)*w0;
eqn = [eqn1 eqn2];
S = solve(eqn,w0,'ReturnConditions',true);
w_0 = double(S.w0);

%% Step 2
r = linspace(R_hub,R_prop,1000);
r_adim = r/R_prop;
phi = atan((V_inf+w_0)./(n_rads*r))';
phi_adim = atan((V_inf+w_0)./(n_rads*r_adim))';
figure
plot(r_adim,phi_adim,'-k');
xlabel('r','interpreter','latex');
ylabel('$\phi$','interpreter','latex');
grid on;

%% Step 3
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

V_eff = (V_inf+w_0)./sin(phi);
figure
plot(r_adim,V_eff,'-k');
xlabel('$\bar{r}$','interpreter','latex');
ylabel('$V_{eff}$','interpreter','latex');
grid on;

lambda = (n_rads*R_prop/V_inf)^-1;
w_j = 2*w_0;
F   = (2/pi)*acos(exp((N_blade/(2*lambda))*(r_adim-1)));
plot(r_adim,F,'-k')
xlabel('$\bar{r}$','interpreter','latex');
ylabel('$F$','interpreter','latex');
grid on;

%Carico
GAMMA = (4*pi.*F.*r_adim.^2.*a_first_radim*n_rads*R_prop^2)/N_blade;
gamma_adim_opti = (N_blade.*GAMMA*n_rads)/(4*pi*V_inf*w_0);
figure
plot(chi,gamma_adim_opti,'-k')
grid on;
xlabel('$\chi$','interpreter','latex');
legend('$\frac{N \Gamma \Omega}{4\pi w_{0}*v_{inf}}$','interpreter','latex');

% 3.31
dT_dr_adim = N_blade*rho*n_rads*(1-a_first_radim).*GAMMA*R_prop^2;
% dP_dr_adim = N_blade*rho*n_rads*(1+a_radim)*V_inf*GAMMA*R_prop^2;

% 3.46
sigma_cl = (1/pi)*((N_blade.*GAMMA)./(V_eff'.*r_adim.*R_prop));

%% Step 4
cl =  [0.8 0.8 0.8 0.6 0.2 0];
r  =  [0.2 0.4 0.6 0.8 0.9 1];
mu = 1e-5;
c = 1;
Re =  (rho.*V_eff*c)/mu;
RE = [Re(200) Re(400) Re(600) Re(800) Re(900) Re(1000)];
alpha = [2.8 2.79 2.77 1 -2.45 -4.19]/57.63;
alpha_int = @(t)interp1(r,alpha,t);
t = linspace(0,1,1000);
alpha_vec = alpha_int(t);
theta = phi(201:1000)+alpha_vec(201:1000);
figure
plot(r_adim(201:1000),theta(:,1),'--k')












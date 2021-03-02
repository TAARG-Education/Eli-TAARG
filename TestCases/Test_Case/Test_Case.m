clear all; close all; clc;
%% Test_case

load("myprop_xrotor.dat");
Cl_xrotor = myprop_xrotor(:,5)';
r_xrotor  = myprop_xrotor(:,2)';

V     = 20;   %m/s
rpm   = 600;
R     = 1.58; %m
r_hub = 0.2;  %m
N     = 3;
ns    = 30;

[alpha_i_deg, beta_deg, c, Cl_alpha] = BEMT(V, rpm, R, r_hub, N, ns);

alpha_i_rad = convang(alpha_i_deg,'deg','rad');
Cl_bemt     = Cl_alpha.*alpha_i_rad;

figure
plot(r_xrotor,Cl_bemt,'k',r_xrotor,Cl_xrotor,'--k');
xlabel('r [m]'); ylabel('C_l');
legend('BEMT','XROTOR');
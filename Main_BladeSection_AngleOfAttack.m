clear all; close all; clc;

%% Data.
Vinf_vec = 1:1:70;
Omega = 220/60*2*pi;
R = 8.50;
A = pi*R^2;
z0 = 2e3;
N = 4;
c = 0.385;
sigma = N*c/(pi*R);
Cl_alpha = 2*pi;
theta_tw = -0.0799;
X = 0;
f_A = 0.015; f = f_A*A;
Lock = 8;
r_segn = linspace(0.15,1,100);
psi = linspace(0,2*pi,100);
alpha_stall_up = 20;
alpha_stall_lo = -10;

%% Calculation.
for i = 1:length(Vinf_vec)
    Vinf = Vinf_vec(i);
    
    [Tc,Hc,Yc,Qc,Pc,alpha_inf,lambda] = Ndim_Coeff_Articulated_Rotor(Vinf,z0,Lock,f,X);
    mu = Vinf*cos(convang(alpha_inf,'deg','rad'))/(Omega*R);
    
    theta0 = 3/(1+3/2*mu^2)*(2*Tc/(sigma*Cl_alpha) - theta_tw/4*(1+mu^2) + lambda/2);
    
    beta0 = Lock*(theta0/8*(1+mu^2)+theta_tw/10*(1+5/6*mu^2)-lambda/6);
    beta1c = -2*mu*(4/3*theta0+theta_tw-lambda)/(1-mu^2/2);
    beta1s = -4/3*mu*beta0/(1+mu^2/2);
    
    beta = beta0 + beta1c*cos(psi) + beta1s*sin(psi);
    dbeta = -beta1c*sin(psi) + beta1s*cos(psi);
    
    theta = theta0 + theta_tw*r_segn;
    
    alpha_stall_up = 20;
    alpha_stall_lo = -10;
    
    alpha_eff = BladeSection_AngleOfAttack(lambda,r_segn,beta,dbeta,mu,psi,theta,alpha_stall_up,alpha_stall_lo);
    
end
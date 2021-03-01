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

function [Tc,Hc,Yc,Qc,Pc,alfa,lambda] = Ndim_Coeff_Articulated_Rotor(V_inf,h,Lock,f,X)
%% Reading geometry and aerodynamics input data
% W,N,D,c,theta_tw,Cl_alfa,Cd_mean and Omega. We think that all these
% parameters should be obtained by a geometry input function, that is not
% yet available. 

%% Input data for local test
% It will be replaced with a function that reads input parameters listed
% below.
 N        = 4;
 W        = 2600;
 R        = 8.50;
 c        = 0.385;
 Cl_alpha = 2*pi;
 Omega    = 220/60*2*pi;
 theta_tw = -0.0799;
 Cd_mean  = 0.01;

%% Data
 [~,~,~,rho_inf] = atmosisa(h);                                             % International Standard Atmosphere model
 g               = 9.8195;                                                  % Gravitational acceleration
 sigma           = N*c/(pi*R);                                              % Rotor solidity
 W               = W*g;                                                      
 A               = pi*R^2;                                                  % Rotor disk area
 D_fus           = f*0.5*rho_inf*V_inf^2;                                   % Fuselage drag force

%% Initialization
 i        = 0;                            % Count iteration
 err      = 1;                            % Error initialization
 err_stop = 1e-5;                         % Error tollerance 
 alfa     = convang(0,'deg','rad');       % Initialization for the angle of attack
 Tc       = W/(rho_inf*(Omega*R)^2*A);    % Thrust coefficient                                      
 X        = convang(X,'deg','rad');       % Rate of climb
 lambda_c = V_inf*sin(X)/(Omega*R);

%% Beginning of the cycle
while abs(err) > err_stop  
    mu = V_inf*cos(alfa)/(Omega*R);                                                           % Rotor advance ratio

%% To overcome the problems due to the divergence of lambda_i at low speed, a different formula provided by Eng. Di Giorgio has been used.    
if mu <= 0.1       
    if i == 0
        lambda_i = sqrt(-V_inf^2/2+sqrt(V_inf^4/4+(W/(2*rho_inf*A))^2))/(Omega*R);            % First attempt value of induced inflow ratio at low speed
    else
        lambda_i = Tc/(2*sqrt(mu^2 + lambda^2));                                              % Updated value of induced inflow ratio 
    end    
end

if mu > 0.1        
    if i == 0
        lambda_i = Tc/(2*mu);                                                                 % First attempt value of induced inflow ratio at high speed
    else
        lambda_i = Tc/(2*sqrt(mu^2 + lambda^2));                                              % Updated value of induced inflow ratio 
    end    
end
%%    
    lambda   = mu*tan(alfa) + lambda_i;                                                     % Rotor inflow ratio
    theta_0  = (2*Tc/(sigma*Cl_alpha) - theta_tw/4*(1 + mu^2) + lambda/2)*3/(1 + 3/2*mu^2); % Collective pitch angle
    Pc0      = sigma*Cd_mean*(1 + 3*mu^2)/8;
    Pc       = lambda_i*Tc + lambda_c*Tc + mu*D_fus*Tc/W + Pc0;                             % Power coefficient

    beta_0   = Lock*(theta_0/8*(1 + mu^2) + theta_tw/10*(1 + 5/6*mu^2) - lambda/6);         % Coning angle
    beta_1c  = -2*mu*(4/3*theta_0+theta_tw-lambda)/(1-mu^2/2);                              % Longitudinal tip-path-plane tilt angle, positive forward
    beta_1s  = -4/3*mu*beta_0/(1+mu^2/2);                                                   % Lateral tip-path-plane tilt angle, positive toward retrating side

    Hci      = sigma*Cl_alpha*0.5*( theta_0*( -1/3*beta_1c+1/2*mu*lambda ) + ...                 % H-induced force coefficient
                 theta_tw*(-1/4*beta_1c + 1/4*mu*lambda) + 3/4*lambda*beta_1c + ...       
                 1/6*beta_0*beta_1s + 1/4*mu*(beta_0^2 + beta_1c^2));
    Hc0      = sigma*Cd_mean*mu/4;                                                               % H-parasite force coefficient   
    Hc       = Hci + Hc0;                                                                        % H-force coefficient
    Qc       = Pc;                                                                               % Torque coefficient
    Yc       = -sigma*Cl_alpha*0.5*( theta_0*(3/4*mu*beta_0 + 1/3*beta_1s*(1 + 3/2*mu^2)) + ...  % Y-force coefficient 
               theta_tw*(1/2*mu*beta_0 + 1/4*beta_1s*(1 + mu^2)) -3/4*lambda*beta_1s      + ...
               beta_0*beta_1c*(1/6 - mu^2)-3/2*mu*lambda*beta_0 - 1/4*beta_1c*beta_1s);

    lambda_old = lambda;                                                                      % Saving old lambda variable                
    lambda     = lambda_i + lambda_c+mu*Hc/Tc + mu*D_fus/W;                                   % Updating lambda at the end of cycle
    alfa       = atan((lambda - Tc/(2*sqrt(mu^2 + lambda^2)))/mu);                            % Updating the value of the angle of attack for the next iteration
    err        = (lambda - lambda_old)/lambda;                                                % Evaluating the error
    i          = i + 1;                                                                       % Increasing the count iteration 
end
%%
alfa = convang(alfa,'rad','deg');
%% Determination of induced and parasite coefficients  
% Qc0 = sigma*Cd_mean*(1 + mu^2)/8; 
% Pc0 = sigma*Cd_mean*(1 + 3*mu^2)/8;
% Qci = Qc - Qc0;
% Pci = Pc - Pc0;
% Yci = Yc;
% Yc0 = 0;
end 
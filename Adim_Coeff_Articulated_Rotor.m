%% Adim_Coeff_Articulated_Rotor.m
%  \brief: the function evaluates the adimensional coefficients of an
%  articulated rotor
%  \author: Marco Artiano, Luca Angelino
%  \version: 1.03
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
% ===============================================================================================
% |Name        : Adim_Coeff_Articulated_Rotor                                                   |
% |Author      : Marco Artiano & Luca Angelino                                                  |
% |              University of Naples Federico II.                                              |
% |Version     : 1.0.3                                                                          |
% |Date        : 25/11/20                                                                       |
% |Modified    : 24/01/21                                                                       |
% |Description : Determination of adimensional aerodynamics coefficients                        |
% |Reference   : 'Lezioni di AERODINAMICA DELL'ALA ROTANTE a.a. 2019-2020 - Renato Tognaccini'  |
% |Input       : Velocity V_inf [m/s], Altitude h [m], Lock number L [~],                       |            
% |              Equivalent drag area f [m^2], Rate of climb X [deg]                            |                                                             
% |Output      : Tc, Hc, Yc, Qc, Pc and the angle of attack [deg]                               |
% |Note        : Missing geometry and aerodynamics input                                        |
% ===============================================================================================

function [Tc,Hc,Yc,Qc,Pc,alfa] = Adim_Coeff_Articulated_Rotor(V_inf,h,Lock,f,X)
%% Reading geometry and aerodynamics input data
% W,N,D,c,theta_tw,Cl_alfa,Cd_mean and Omega. We think that all these
% parameters should be obtained by a geometry input function, that is not
% yet available. 

%% Input data for local test
% It will be replaced with a function that reads input parameters listed
% below.
 N        = 4;
 W        = 2600;
 R        = 11.9/2;
 c        = 0.385;
 Cl_alpha = 2*pi;
 Omega    = 40.3361;
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
if V_inf <= 20       
    if i == 0
        lambda_i = sqrt(-V_inf^2/2+sqrt(V_inf^4/4+(W/(2*rho_inf*A))^2))/(Omega*R);            % First attempt value of induced inflow ratio at low speed
    else
        lambda_i = Tc/(2*sqrt(mu^2 + lambda^2));                                              % Updated value of induced inflow ratio 
    end    
end
       
if V_inf > 20        
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
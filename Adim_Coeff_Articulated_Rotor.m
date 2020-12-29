function [Tc,Hc,Yc,Qc,Pc,alfa] = Adim_Coeff_Articulated_Rotor(V_inf,altitude,Lock,f,X)
%% Adim_Coeff_Articulated_Rotor
% Determination of adimensional aerodynamics coefficients: 
% - Tc   thrust coefficient 
% - Hc   rotor drag force coefficient 
% - Yc   lateral force coefficient 
% - Qc   torque coefficient
% - Pc   power coefficient
% - alfa angle of attack
%
% Marco Artiano - Luca Angelino
% version 1.0.0
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
% |Version     : 1.0.0                                                                          |
% |Date        : 25/11/20                                                                       |
% |Modified    : 06/12/20                                                                       |
% |Description : Determination of adimensional aerodynamics coefficients                        |
% |Reference   : 'Lezioni di AERODINAMICA DELL'ALA ROTANTE a.a. 2019-2020 - Renato Tognaccini'  |
% |Input       : Velocity with respect to the air (V_inf), Altitude(altitude)                   |
% |Output      : Tc, Hc, Yc, Qc, Pc and CM                                                      |
% |Note        : Missing geometry and aerodynamics input and determination of moment coefficient|
% ===============================================================================================

%% Input parameters
%  V_inf     [m/s]     Rotor or helicopter velocity with respect to the air
%  altitude  [m]       Altitude
%  W         [kg]      Helicopter or rotor gross weight
%  N         [~]       Number of blades
%  D         [m]       Rotor diameter
%  c         [m]       Rectangular blade chord
%  theta_tw  [rad]     Linear twist rate
%  Cl_alpha  [1/rad]   Lift coefficient slope
%  Cd_mean   [~]       Mean drag coefficient 
%  Omega     [rad/s]   Rotor rotational speed
%  f         [m^2]     Equivalent drag area of helicopter fuselage and hub
%  Lock      [~]       Lock number
%  X         [deg]     Rate of climb 
 
%% Note
% W,N,D,c,theta_tw,Cl_alfa,Cd_mean,Omega,f and Lock. We think that all these
% parameters should be obtained by a geometry input function, that is not
% yet available. 

%% Input data for local test
% It will be replaced with a function that reads input parameters listed below
 N        = 4;
 W        = 2600;
 R        = 11.9/2;
 c        = 0.385;
 Cl_alpha = 2*pi;
 Omega    = 40.3361;
 theta_tw = -0.0799;
 Cd_mean  = 0.01;
 A        = pi*R^2; 
 
 %% Test values  
 if nargin == 0
 f        = 0.007*A;
 Lock     = 8;
 X        = convang(18,'deg','rad'); 
 V_inf    = 20;
 altitude = 0;
 end

%% Reading geometry and aerodynamics input data
%  not yet available 
%
%% Data
 [~,~,~,rho_inf] = atmosisa(altitude);                                      % International Standard Atmosphere model
 g               = (6.674*10^(-11))*(5.972*10^24)/(6371000+altitude)^2;     % Gravitational acceleration with respect to the altitude
 sigma           = N*c/(pi*R);                                              % Rotor solidity
 W               = W*g;                                                      
 A               = pi*R^2;                                                  % Rotor disk area
 D_fus           = f*0.5*rho_inf*V_inf^2;                                   % Fuselage drag force
 
%% Initialization
 i        = 0;                            % Count iteration
 err      = 1;                            % Error initialization
 err_stop = 1e-3;                         % Error tollerance 
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

%% Determination of induced and parasite coefficients  
% Qc0 = sigma*Cd_mean*(1 + mu^2)/8; 
% Pc0 = sigma*Cd_mean*(1 + 3*mu^2)/8;
% Qci = Qc - Qc0;
% Pci = Pc - Pc0;
% Yci = Yc;
% Yc0 = 0;

%% Plotting inputs and outputs
 table_disp_in = table([V_inf; altitude; W/g; N; 2*R; c; theta_tw; Cl_alpha; Cd_mean; Omega; f/A; Lock],...
                 categorical(["m/s";" m ";" kg ";" ~ ";" m ";" m ";" rad ";" 1/rad ";" ~ ";" rad/s ";" ~ ";" ~ "]),...
                 'VariableNames',{'Value','Measure unit'},...
                 'RowName',{'V_inf';'Altitude';'W';'N';'D';'c';'theta_tw';'Cl_alpha';'Cd_mean';'Omega';'f/A';'Lock'});
            
 table_disp_ou = table([Tc; Hc; Yc; Qc; Pc; alfa*57.3],...
                 'VariableNames',{'Value'},...
                 'RowName',{'Tc';'Hc';'Yc';'Qc';'Pc';'alfa [deg]'});
 fprintf('\n');
 disp('<strong>Input data:</strong>');
 disp(table_disp_in);
 disp('<strong>Output data:</strong>');
 disp(table_disp_ou);

end
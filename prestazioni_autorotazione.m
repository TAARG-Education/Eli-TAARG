%% \file prestazioni_autorotazione
%  \brief Autorotation performances
%  \author Giulio Guarino - Emanuele Soreca
%  \version 1.0.0
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
% |Name        : prestazioni_autorotazione                                                     |
% |Author      : Giulio Guarino - Emanuele Soreca                                              |
% |              University of Naples Federico II.                                             |
% |Version     : 1.0.0                                                                         |
% |Date        : 13/01/2020                                                                    |
% |Modified    : 18/01/2020                                                                    |
% |Description : This function provides autorotation performances                              |
%                of a rotor: the procedure is based upon the one shown                         | 
%                in reference. For a given helicopter, advance ratio                           | 
%                and descent angle are to be assigned by the user.                             | 
%                Autorotation condition yields to equations (see                               |
%                reference, page 110) which are here dealt with                                | 
%                by means of the Matlab Symbolic Toolbox. Once solved                          | 
%                in respect to the intake ratio, expressions of thrust                         |
%                coefficient, angle of attack, angular velocity, tip                           | 
%                velocity, freestream velocity can be evaluated.                               | 
% |Reference   : Lezioni di Aerodinamica dell'Ala Rotante a.a. 2019/2020                       |
%                Renato Tognaccini                                                             | 
% |Input       : Advance ratio, descent ratio values (to be assigned).                         |
% |Output      : Numerical values of thrust coefficient, angle of                              |
%                attack, angular velocity, freestream velocity, tip                            |
%                velocity.                                                                     | 
% |Note        : note addizionali                                                              |
% ==============================================================================================

function [T_C, H_C, alpha_deg, omega, V_inf ] = prestazioni_autorotazione(mu,X)

global  R N c gamma theta_tw f rho M

Cl_a = 2*pi;                          % [1/rad], lift coefficient gradient (from thin airfoil theory)
Cd = 0.01;                            % Average drag coefficient along the blade
X = deg2rad(X);                       % [rad], descent angle
lambda_c = -mu*sin(X);                % Descent ratio (<0)
A = 3.14*R^2;                         % [m^2], swept area
sigma = N*c/(pi*R);                   % rotor solidity
W = M*9.81;                           % [N], weight


%% Beginning of the procedure.
% We need to evaluate an expression of theta_0(lambda). Evaluate T_c, Q_c0,
% H_ci first:

syms theta_0 lambda real

T_c = 0.5*sigma*Cl_a*(theta_0/3*(1+3/2*mu^2)+theta_tw/4*(1+mu^2)-lambda/2);  % thrust coefficient

Q_c0 = sigma*Cd/8*(1+mu^2);                                                  % parasitic torque coefficient

b_0 = gamma*(theta_0/8*(1+mu^2)+theta_tw/10*(1+5/6*mu^2)-lambda/6);          % flapping coefficients expressions

b_1c = -2*mu*(4/3*theta_0+theta_tw-lambda)/(1-0.5*mu^2);

b_1s = -4/3*mu*b_0/(1+0.5*mu^2);

H_ci = 0.5*sigma*Cl_a*(theta_0*(-1/3*b_1c+0.5*mu*lambda)...                  % induced rotor drag coefficient
    +theta_tw*(-1/4*b_1c+1/4*mu*lambda)+3/4*lambda*b_1c...
    +1/6*b_0*b_1s+1/4*mu*(b_0^2+b_1c^2));

% We can find an expression of theta_0(lambda) then (we use capital
%subscripts from this point)

%from the definition of autorotation:
P_c = lambda*T_c+Q_c0-mu*H_ci;                                               % Equation to solve in respect to theta_0

THETA_0 = solve(P_c==0,theta_0);                                             % Two solutions

THETA_0 = THETA_0(2);                                                        % We choose the second one

% It is possible now to find Tc=Tc(lambda),Hc=Hc(lambda),alfa=alfa(lambda)

T_C = 0.5*sigma*Cl_a*(THETA_0/3*(1+3/2*mu^2)...                             % Tc as a function of lambda only
    +theta_tw/4*(1+mu^2) -lambda/2); 

H_C0 = sigma*Cd*mu/4;                                                       % parasitic rotor drag coefficient

B_0 = gamma*(THETA_0/8*(1+mu^2)+theta_tw/10*(1+5/6*mu^2)...                 % flapping coefficients as a function of lambda only
    -lambda/6);    

B_1c = -2*mu*(4/3*THETA_0+theta_tw-lambda)/(1-0.5*mu^2);

B_1s_lambda = -4/3*mu*B_0/(1+0.5*mu^2);

H_Ci = 0.5*sigma*Cl_a*(THETA_0*(-1/3*B_1c...                                 % induced rotor drag coefficient
    +0.5*mu*lambda)+theta_tw*(-1/4*B_1c+1/4*mu*lambda)...                    % as a function of lambda only
    +3/4*lambda*B_1c+1/6*B_0*B_1s_lambda...
    +1/4*mu*(B_0^2+B_1c^2));

alpha = atan(1/mu*(lambda-T_C/(2*sqrt(mu^2+lambda^2))));                    % angle of attack as a function of lambda only

% We can now solve an equation for lambda

lambda_i = T_C/(2*mu);
%lambda_i = T_C/(2*sqrt(mu^2+lambda^2))

P_C = T_C*lambda_i+lambda_c*T_C+...                                         % This expression is obtained from the Pc expression
      mu*0.5*f/A*(mu/cos(alpha))^3+H_C0+Q_c0;

lambda = vpasolve(P_C==0,lambda);

%% Evaluating output

T_C       = double(subs(T_C))                                                   % Thrust coefficient

H_C       = double(subs(H_Ci))+H_C0                                             % Rotor drag coefficient

alpha     = double(subs(alpha));                                                % [rad] Angle of attack

omega     = double(sqrt(W/(rho*A*T_C*R^2)))                                     % [rad/s] Angular velocity

V_inf     = mu*omega*R/cos(alpha)                                               % [m/s] Freestream velocity

alpha_deg = convang(alpha,'rad','deg')                                          % [deg] Angle of attack in degrees
end



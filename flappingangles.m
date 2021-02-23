%% \flappingangles.m
%  \brief: the function evaluates the flapping angles of a rotor blade
%  in forward flight
%  \authors: Matteo Cassetta, Pasquale Antonio Federico
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
% |Name        : flappingangles.m                                                        |
% |Authors     : Matteo Cassetta                                                         |
%                Pasquale Antonio Federico                                               |
% |              University of Naples Federico II                                        |
% |Version     : 1.0                                                                     |
% |Date        : 08/01/2021                                                              |
% |Modified    : 26/01/2021                                                              |
% |Description : Rotor blade flapping angles                                             |
% |Reference   : Prouty, R.W. - Helicopter Performance, Stability and Control            |
% |Input       : Airspeed, angle of attack (deg), rotor/helicopter gross weight,         |
% |              rotor Lock number, rotor solidity,rotor blade twist (deg,assumed linear)|
% |              blade airfoil lift curve slope (1/rad) , rotor blade tip speed,         |
% |              flapping hinge eccentricity/blade radius, blade radius;                 |
% |              optional inputs allow choosing of the output mode                       |
% |              and azimuth angle sample in "angles" mode                               |
% |Output      : Azimuth angles psi, flapping angles beta and dbeta/dpsi                 |
% |Note        :                                                                         |
% ==============================================================================================

function [varargout] = flappingangles(V,alpha,W,gamma,...
    sigma,theta_tw,cla,omegaR,eR,R,options)

% Argument validation
arguments
    V        (1,1) {mustBeNumeric,mustBeReal}
    alpha    (1,1) {mustBeNumeric,mustBeReal}
    W        (1,1) {mustBeNumeric,mustBeNonnegative}
    gamma    (1,1) {mustBeNumeric,mustBeNonnegative}
    sigma    (1,1) {mustBeNumeric,mustBeNonnegative}
    theta_tw (1,1) {mustBeNumeric,mustBeReal}
    cla      (1,1) {mustBeNumeric,mustBeReal}
    omegaR   (1,1) {mustBeNumeric,mustBeNonnegative}
    eR       (1,1) {mustBeNumeric,mustBeNonnegative}
    R        (1,1) {mustBeNumeric,mustBeNonnegative}
    options.output  (1,1) string {mustBeMember(options.output,{'angles','coefficients'})} = 'angles'
    options.sample  (:,1) {mustBeNumeric,mustBeReal} = linspace(0,360,200)'
    
end

%% Preliminary calculations
altitude = 0;
[~,~,~,rho] = atmosisa(altitude);

% Converting to radians
theta_tw = convang(theta_tw,'deg','rad');

% Thrust coefficient
Tc = W/(rho*omegaR^2*pi*R^2);

% Advance ratio
mi = V*cosd(alpha)/omegaR;

% Adimensional inflow
lambda_i = sqrt(-.5*V.^2+.5*sqrt(V.^4+4*(W/(2*rho*pi*R^2)).^2))/omegaR;
lambda   = mi.*tand(alpha)+ lambda_i;

% Collective pitch
theta_0 = 6*Tc./(sigma*cla*(1+3/2*mi.^2)*(1-eR))...
    -3/4*theta_tw*(1+mi.^2)./(1+3/2*mi.^2)+1.5*lambda./(1+3/2*mi.^2);

%% Flapping coefficients

% Coning
beta_0 = 1/6*gamma*(3/4*theta_0.*(1+mi.^2)+theta_tw*(3/5+mi.^2/2)-lambda)...
    *(1-eR)^2;

% Longitudinal flapping
beta_1c = -2*mi.*(4/3*theta_0+theta_tw-lambda)./(1-mi.^2/2)...
    -12*eR./(gamma*(1-eR)^3*(1+mi.^4/4)).*(4/3*Tc/sigma*...
    (2/3*mi*gamma/(cla*(1+3/2*eR))));

% Lateral flapping
beta_1s = -4/3*mi.*beta_0.*(1+eR/2)./((1+mi.^2/2)*(1-eR)^2)...
    -12*eR*(1+eR/2)./(gamma*(1-eR)^3*(1+mi.^4/4))*...
    2.*mi.*(4/3*theta_0+theta_tw-lambda);

%% Output

if strcmpi(options.output,'coefficients')
    varargout = {convang(beta_0,'rad','deg'),...
        convang(beta_1c,'rad','deg'),...
        convang(beta_1s,'rad','deg')};
else
    psi    = options.sample;
    beta = (beta_0 + beta_1c.*cosd(psi) + beta_1s.*sind(psi));
    beta_dot = ( -beta_1c.*sind(psi) + beta_1s.*cosd(psi));
    varargout = {psi,convang(beta,'rad','deg'),...
        convang(beta_dot,'rad','deg')};
end

end
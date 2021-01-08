%% \flappingangles.m
%  \brief: the function evaluates the flapping angles of a rotor blade
%  in forward flight
%  \authors: Matteo Cassetta, Pasquale Antonio Federico
%  \version: 1.04
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
% |Modified    : 08/01/2021                                                              |
% |Description : Rotor blade flapping angles                                             |     
% |Reference   : Prouty, R.W. - Helicopter Performance, Stability and Control            |
% |Input       : Airspeed, angle of attack (rad), rotor/helicopter gross weight,         |
% |              rotor Lock number, rotor solidity,rotor blade twist (rad,assumed linear)|
% |              blade airfoil lift curve slope (1/rad) , rotor blade tip speed,         |  
% |              flapping hinge eccentricity/blade radius, blade radius;                 |
% |              optional inputs allow choosing of the output mode                       |
% |              and azimuth angle sample in "angles" mode                               |
% |Output      : Azimuth angles psi, flapping angles beta and dbeta/dpsi                 |
% |Note        :                                                                         |
% ==============================================================================================

function [varargout] = flappingangles(V,alpha,W,gamma,...
                            sigma,theta_tw,cla,omegaR,eR,R,varargin)

validInputs = {'output','sample'};

% Default sample and output
output = 'angles';
psi = linspace(0,2*pi,200);

%% Optional input check

for i = 1:length(varargin)
    if mod(i,2) == 1
        % Cast uneven index inputs to string to check their validity as
        % keys
        varargin{i} = num2str(varargin{i});
    end
    % Converting to lowercase to make the control case insensitive
    varargin{i} = lower(varargin{i});
end

% A total of 2 optional input arguments + values are allowed.
% Any further input arguments are discarded.
if length(varargin) > 4
    warning('Too many optional inputs, some will be ignored.')
    varargin{5:end} = [];
        
else
    % If there is an uneven number of input it is rendered even by removing
    % the last input.
    if  mod(length(varargin),2) == 1
        if any(ismember(varargin(end),validInputs))
            warning(['No input argument found for "',varargin{end},...
                '".']);   
        else
            warning([varargin{end},' is not a valid optional input and will be ignored.']);
        end
        varargin(end) = [];
    end

    % Optional inputs are divided in pairs: input,argument
    for i = 1:2:length(varargin)
            % Check if input key is valid
            if ismember(varargin(i),validInputs)
                % Check if argument is valid then assign it
                if strcmpi(varargin(i),'output')
                    if any(strcmpi(varargin{i+1},{'angles','coefficients'}))
                        output = varargin{i+1};
                    else
                        warning([num2str(varargin{i+1}), ' is not a valid argument for "',...
                            varargin{i},'" and will be ignored.']);
                    end
                else
                    if isnumeric(varargin{i+1})
                        psi = varargin{i+1};
                    else
                        warning(['Non numerical input for "',varargin{i},...
                            '" is not valid and will be ignored.']);
                    end
                end
            else
                warning([varargin{i},' is not a valid optional input and will be ignored.']);
            end
    end
end

% Thrust coefficient
Tc = W/(1.225*omegaR^2*pi*R^2);

% Advance ratio
mi = V*cos(alpha)/omegaR;

% Adimensional inflow
lambda_i = sqrt(-.5*V.^2+.5*sqrt(V.^4+4*(W/(2*1.225*pi*R^2)).^2))/omegaR;
lambda   = mi.*tan(alpha)+ lambda_i;

% Collective pitch
theta_0 = 6*Tc./(sigma*cla*(1+3/2*mi.^2)*(1-eR))...
    -3/4*theta_tw*(1+mi.^2)./(1+3/2*mi.^2)+1.5*lambda./(1+3/2*mi.^2);


%%% Flapping angle coefficients

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
        
if strcmpi(output,'coefficients')
    varargout = {beta_0,beta_1c,beta_1s};
else
    beta = beta_0 + beta_1c*cos(psi) + beta_1s*sin(psi);
    beta_dot = beta_1c*sin(psi) - beta_1s*cos(psi);
    varargout = {psi,beta,beta_dot};
end

end
%% \file PresRot.m
%  \brief: Rotor performance analysis in hoverig and climbing conditions
%  \author: Gianpietro Cinque, Verdiana Bottino
%  \version: 1.00
%
% Eli-TAARG is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public
% License as published by the Free Software Foundation; either
% version 3 of the License, or (at your option) any later version.
%
% Eli-TAARG is developed by the TAARG Educational organization for
% educational purposes only.
% Theoretical and Applied Aerodynamic Research Group - University of
% Naples Federico II.
%
% Eli-TAARG GitHub link: <https://github.com/TAARG-Education/Eli-TAARG>
%
% Eli-TAARG is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
% General Public License for more details.
% <http://www.gnu.org/licenses/>.
%
% ========================================================================================
% |Name        : PresRot.m                                                               |
% |Author      : Verdiana Bottino, Gianpietro Cinque                                     |
% |              University of Naples Federico II.                                       |
% |Version     : 1.00                                                                    |
% |Date        : 24/04/2021                                                              |
% |Modified    : 29/04/2021                                                              |
% |Description : Given the main geometric and aerodynamic parameters, the
% |function retrives thrust and torque coefficients using the blade element
% |theory
% |Reference   : "Teoria del volo dell'elicottero" - Giovanni di Giorgio,
% |              Aracne Editrice
% |                                                                                      |                                                           |
% |Input       : Rr [m] - hub radius
% |              Rt [m] - tip radius
% |              thetaR [°] - collective pitch
% |              deltaTheta [°] - rotor geometric twist
% |              omega [rad/s] - angular velocity
% |              chord [m]
% |              N - number of blades
% |              Vc [m/s] - rate of climb
% |              rho [kg/m^3]- air density
% |              airfoil - airfoil 'filename'; this file must be placed in
% |              the running directory and consistent with Xfoil input
% |              requirements                                                            |
% |                                                                                      |
% |                                                                                      |
% |                                                                                      |
% |                                                                                      |
% |                                                                                      |
% |                                                                                      |
% |                                                                                      |
% |                                                                                      |
% |                                                                                      |
% |                                                                                      |
% |Output      : Tc - thrust coefficient
% |              Qc - torque coefficient
% |                                                                                      |
% |                                                                                      |
% |Note        :                                                                         |
% ========================================================================================

function [ Tc, Qc ] = PresRot( Rr, Rt, thetaT, deltaTheta,omega, chord, N, Vc, rho, airfoil, varargin )


%--------------------------------------------------------------------------
% VARIABLES
%--------------------------------------------------------------------------
Ns = 100;                                        % Number of blades
thetaVet = linspace( thetaT, thetaT - deltaTheta, Ns );
rVet = linspace( Rr / Rt, 1, Ns );
mu = Vc / ( omega * Rt );
Cla = 2 * pi;
Np = 200;                                        % Number of panels, Xfoil
viscosity = 1.78e-5;

%--------------------------------------------------------------------------
% VARIABLE INPUT
%--------------------------------------------------------------------------
n = length( varargin );

if isempty( varargin )
    cVet = chord * ones( 1, Ns );
else
    
end

%--------------------------------------------------------------------------
% ALLOCATE VARIABLES
%--------------------------------------------------------------------------
sigmaVet = N * cVet / ( 2 * pi * Rt );
lambdaVet=zeros(1,Ns);
phiVet=zeros(1,Ns);
alphaVet=zeros(1,Ns);
ReVet = zeros( 1, Ns );
ClVet=zeros(1,Ns);
CdVet=zeros(1,Ns);
dTdr=zeros(1,Ns);
dQdr=zeros(1,Ns);

%--------------------------------------------------------------------------
% Checks
%--------------------------------------------------------------------------
if ( not( exist( 'rho', 'var' ) ) )
    rho = 1.225;
end

%--------------------------------------------------------------------------
% Calculate Tc & Qc
%--------------------------------------------------------------------------

for i = 1 : Ns
    
    theta = thetaVet( i );
    r = rVet( i );
    chordi = cVet( i );
    sigma = sigmaVet( i );
    
    a = 1;
    b = mu + Cla *(sigma/8);
    c = -(r * Cla *(sigma/8)*(theta-(mu/r)));
    
    % Assemble lambda vector, accept only the positive root
    lambdaVet(i) = max( roots( [ a b c ] ) );
    phiVet(i) = lambdaVet(i)/ rVet(i);         % inflow angle
    alpha = theta - phiVet(i);                 % effective aoa
    alphaVet(i) = alpha;
    
    % Effective velocity computation through the blade element theory
    Ve = sqrt( Vc^2 + ( omega * r * Rt )^2 );
    Ve = ( omega * r ) / ( cosd( phiVet( i ) ) );
    % Reynolds number computation
    Re = round( ( rho * Ve * chordi ) / ( viscosity ) );
    ReVet( i ) = Re;
    
    [ ClVet( i ), CdVet( i ) ] = ...
        CdCl_xfoil( airfoil, Np, Re, alpha, alpha, 1 );
    
    dTdr(i)=1/2*sigma*ClVet(i)*rVet(i)^2;
    dQdr(i)=1/2*sigma*(ClVet(i)*phiVet(i)+CdVet(i))*(rVet(i)^3);
    
    
end

Tc = sum( dTdr );
Qc = sum( dQdr );

end

function [ Ni ] = inputN ( name )

switch name
    case 'variableChord'
        Ni = 2;
    case 'variableTheta'
        Ni = 2;
    otherwise
        error(message('MATLAB:PresRot:wrongVariableName'));
end

end






%--------------------------------------------------------------------------
% Emprovements
%--------------------------------------------------------------------------

% si potrebbe calcolare il Cla per i lambdai tramite Xfoil
% variazione di rho con la quota
% variazione della viscosty in funzione della quota
% miglior approssimazione della Veffettiva
% in CdCl_xfoil manca il comando per uscire alla fine di Xfoil
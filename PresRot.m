%% \file PresRot.m
%  \brief: Analysis of the rotor performances in hoverig and ??? conditions
%  \author: Verdiana Bottino, Gianpietro Cinque
%  \version: 1.00
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
% |Name        : PresRot.m                                                           |
% |Author      : Verdiana Bottino, Gianpietro Cinque                                                        |
% |              University of Naples Federico II.                                       |
% |Version     : 1.00                                                                    |
% |Date        : 24/04/2021                                                              |
% |Modified    : 24/04/2021                                                              |
% |Description : Given the ......., the function derives and retrives the.....           |     
% |Reference   : "Teoria del volo dell'elicottero" - Giovanni di Giorgio, 
% |              Aracne Editrice
% |                          |                                                         |
% |Input  ????     :                                                                         |
% |                                                                                  |
% |                              |  
% |                                           |
% |                                            |
% |                                                 |
% |                                                |
% |                                   |
% |                                   |
% |                | 
% |                                |
% |Output      :                                                                         |
% |                                                                                      |
% |                                                                                      |
% |Note        :                                                                         |
% ==============================================================================================

function [ Tc, Qc ] = PresRot( Rr, Rt, thetaT, deltaTheta,omega, chord, N, Vc, rho, airfoil, vargin )
% eventualmente: corda variabile, theta non lineare

%--------------------------------------------------------------------------
% Variables
%--------------------------------------------------------------------------
Ns = 100;                                    % Number of stations
thetaVet = linspace( thetaT, thetaT - deltaTheta, Ns );%
rVet = linspace( Rr / Rt, 1, Ns );
mu = Vc / ( omega * Rt );
Cla = 2 * pi;
sigma = N * chord / ( 2 * pi * Rt );
Np = 200;                                    % Number of panels
viscosity = 1.78e-5;                                 % implementa
%--------------------------------------------------------------------------
% Allocate variables
%--------------------------------------------------------------------------
lambdaVet=zeros(1,Ns);
phiVet=zeros(1,Ns);
alphaVet=zeros(1,Ns);
ReVet = zeros( 1, Ns );
ClVet=zeros(1,Ns);
CdVet=zeros(1,Ns);
TcVet=zeros(1,Ns);
QcVet=zeros(1,Ns);
dTdr=zeros(1,Ns);
dQdr=zeros(1,Ns);

%--------------------------------------------------------------------------
% Check
%--------------------------------------------------------------------------
if ( not( exist( 'rho', 'var' ) ) )
    rho = 1.225;
end

%--------------------------------------------------------------------------
% Find Tc & Qc
%--------------------------------------------------------------------------

for i = 1 : Ns
    
    theta = thetaVet( i );
    r = rVet( i );
    
    a = 1;
    b = mu + Cla *(sigma/8);
    c = -(r * Cla *(sigma/8)*(theta-(mu/r)));
    
    % Assemble lambda vector
    lambdaVet(i) = max( roots( [ a b c ] ) );% accept only the positive root
    phiVet(i) = lambdaVet(i)/ rVet(i);
    alpha = theta - phiVet(i);
    alphaVet(i) = alpha;
    
    % Calcolo di un Re approssimato utilizzanto la teoria dell'elemento di
    % pala semplice
    Ve = sqrt( Vc^2 + ( omega * r * Rt )^2 );
    Ve = ( omega * r ) / ( cosd( phiVet( i ) ) );
    Re = ( rho * Ve * chord ) / ( viscosity );
    ReVet( i ) = Re;
    
    [ ClVet( i ), CdVet( i ) ] = CdCl_xfoil( 4412, Np, Re, alpha, alpha, 1 );
    
    TcVet(i)=2*lambdaVet(i)^2;
    QcVet(i)=(TcVet(i)^(3/2))/sqrt(2);
    dTdr(i)=1/2*sigma*ClVet(i)*rVet(i)^2;
    dQdr(i)=1/2*sigma*(ClVet(i)*phiVet(i)+CdVet(i))*(rVet(i)^3);
    
    
end
plot( TcVet, QcVet, 'k' );

Tc = sum( dTdr );
Qc = sum( dQdr );

%--------------------------------------------------------------------------
% Emprovements
%--------------------------------------------------------------------------

% si potrebbe calcolare il Cla per i lambdai tramite xfoil
% numeri di Mach e correzione Glauert?
% rho che varia con la quota
% viscosty in funzione della quota
% miglior approssimazione della Veffettiva
end


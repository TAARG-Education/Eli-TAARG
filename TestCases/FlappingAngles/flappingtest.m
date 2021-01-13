%% \flapping_test.m
%  \brief: this script test the flappingangles.m function
%  \authors: Matteo Cassetta, Pasquale Antonio Federico
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
% |Name        : flapping_test.m                                                         |
% |Authors     : Matteo Cassetta                                                         |
%                Pasquale Antonio Federico                                               |
% |              University of Naples Federico II                                        |
% |Version     : 1.0                                                                     |
% |Date        : 08/01/2021                                                              |
% |Modified    : 13/01/2021                                                              |
% |Description : Rotor blade flapping angles                                             |     
% |Reference   : Prouty, R.W. - Helicopter Performance, Stability and Control            |
%                Johnson, W.  - Helicopter Theory                                        |
% |Input       : Airspeed, angle of attack (rad), rotor/helicopter gross weight,         |
% |              rotor Lock number, rotor solidity,rotor blade twist (rad,assumed linear)|
% |              blade airfoil lift curve slope (1/rad) , rotor blade tip speed,         |  
% |              flapping hinge eccentricity/blade radius, blade radius;                 |
% |              optional inputs allow choosing of the output mode                       |
% |              and azimuth angle sample in "angles" mode                               |
% |Output      : Flapping angle coefficient curves                                       |
% |Note        :                                                                         |
% ==============================================================================================

clear all; close all; clc;

% Declare velocity array and flapping coefficients arrays
[b0,b1c,b1s,v] = deal(linspace(0,120,200));

% Calculation
for i = 1:length(v)
    [b0(i),b1c(i),b1s(i)] = flappingangles(v(i),0,33523.27,8,0.1,-8,5.7,213,0,4,...
        'output','coefficients');
end

% Plots
figure; subplot(2,2,[1 2]);
plot(v,b0,'-k');
xlim([v(1) v(end)])
xlabel('$V$','Interpreter','latex','fontsize',14);
ylabel('$\beta_0 \ [\circ]$','Interpreter','latex','fontsize',14);

subplot(2,2,3);
plot(v,-b1c,'-k');
xlim([v(1) v(end)])
xlabel('$V$','Interpreter','latex','fontsize',14);
ylabel('$-\beta_{1c} \ [\circ]$','Interpreter','latex','fontsize',14);
 
subplot(2,2,4);
plot(v,-b1s,'-k');
xlim([v(1) v(end)])
xlabel('$V$','Interpreter','latex','fontsize',14);
ylabel('$-\beta_{1s} \ [\circ]$','Interpreter','latex','fontsize',14);




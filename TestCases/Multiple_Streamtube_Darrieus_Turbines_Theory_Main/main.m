%% \main.m
%  \brief: this script tests the Multiple_Streamtube_Darrieus_Turbines_Theory.m function
%  \author: Luciano Roncioni - Camilla Scotto di Carlo
%  \version: 1.1
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
% |Name        : main.m                                                                        |
% |Author      : Luciano Roncioni - Camilla Scotto di Carlo                                    |
% |              University of Naples Federico II.                                             |
% |Version     : 1.1                                                                           |
% |Date        : 28/12/2020                                                                    |
% |Modified    : 28/12/2020                                                                    |
% |Description : turbine_Darrieus_tubo_flusso_multiplo.m testing                               |     
% |Reference   : Tognaccini R., (2019), "Lezioni di Aerodinamica dell'ala rotante"             |
%                De Vries O., (1979), "Fluid Dynamic Aspects of Wind Energy Conversion"        |                                                       
% |Input       : A vector with (alphamax,c,R,cla,N,cd) for four                                |
%                  different cases, where :                                                    |
%                 (alphamax) = angle of attack at the stall of the considered airfoil          |
%                 (c) = chord of the blade                                                     |
%                 (R) = blade radius                                                           |
%                 (Cla) = slope coefficient of the Cl curve for the considered airfoil         |
%                 (N) = number of blades                                                       |
%                 (Cd) = drag coefficient of the considered airfoil                            |
%                 (DeVries_data.mat) = a .mat file containing the                              |
%                 corresponding test case                                                      |
%                 literature values from De Vries, "Fluid Dynamic Aspects                      | 
%                 of Wind Energy Conversion"                                                   |
% |Output      :  (cp) = vector of power coefficient values for the                            |
%                 turbine                                                                      |
%                 (cq) = vector of torque coefficient values for the                           |
%                 turbine                                                                      |
%                 (lambdav) = vector of tip speed values corresponding to                      |
%                 the cp and cq coefficients                                                   |
% |Note        :  results have been validated through literature values                        |
% ==============================================================================================
clear all; clc; close all;

% This is the main code written to validate the function.

%Test Cases values
test_n1 = [14,1,30,6.28,3,0.01];
test_n2 = [14,1,30,6.28,3,0];
test_n3 = [14,1,15,6.28,3,0.01];
test_n4 = [14,1,15,6.28,3,0];

load DeVries_data.mat %import data from De Vries

%% Cases
n = 1;
switch n
    case 1
        alphamax = test_n1(1);
        c = test_n1(2);
        R = test_n1(3);
        cla = test_n1(4);
        N = test_n1(5);
        cd = test_n1(6);
        %% Function
        [cp,cq,lambdav] = Multiple_Streamtube_Darrieus_Turbines_Theory(alphamax,c,R,cla,N,cd);
        %% Plot
        plot(lambda1,cp1,'kv'); hold on; grid on;
        plot(lambdav,cp,'k-');  
        legend('single streamtube theory','multiple streamtube theory','Location','best');
        h=xlabel('$\lambda$','Interpreter','latex'); set(h,'Fontsize',20);
        h=ylabel('$C_P$','Interpreter','latex');
        set(h,'Fontsize',20);
        title('\sigma = 0.1 - Cd = 0.01');
        
    case 2
        alphamax = test_n2(1);
        c = test_n2(2);
        R = test_n2(3);
        cla = test_n2(4);
        N = test_n2(5);
        cd = test_n2(6);
        %% Function
        [cp,cq,lambdav] = Multiple_Streamtube_Darrieus_Turbines_Theory(alphamax,c,R,cla,N,cd);
        %% Plot
        plot(lambda2,cp2,'kv'); hold on; grid on;
        plot(lambdav,cp,'k-');  
        legend('single streamtube theory','multiple streamtube theory','Location','best');
        h=xlabel('$\lambda$','Interpreter','latex'); set(h,'Fontsize',20);
        h=ylabel('$C_P$','Interpreter','latex');
        set(h,'Fontsize',20);
        title('\sigma = 0.1 - Cd = 0');
        
    case 3
        alphamax = test_n3(1);
        c = test_n3(2);
        R = test_n3(3);
        cla = test_n3(4);
        N = test_n3(5);
        cd = test_n3(6);
        %% Function
        [cp,cq,lambdav] = Multiple_Streamtube_Darrieus_Turbines_Theory(alphamax,c,R,cla,N,cd);
        %% Plot
        plot(lambda3,cp3,'kv'); hold on; grid on;
        plot(lambdav,cp,'k-');  
        legend('single streamtube theory','multiple streamtube theory','Location','best');
        h=xlabel('$\lambda$','Interpreter','latex'); set(h,'Fontsize',20);
        h=ylabel('$C_P$','Interpreter','latex');
        set(h,'Fontsize',20);
        title('\sigma = 0.2 - Cd = 0.01');
        
    case 4
        alphamax = test_n4(1);
        c = test_n4(2);
        R = test_n4(3);
        cla = test_n4(4);
        N = test_n4(5);
        cd = test_n4(6);
        %% Function
        [cp,cq,lambdav] = Multiple_Streamtube_Darrieus_Turbines_Theory(alphamax,c,R,cla,N,cd);
        %% Plot
        plot(lambda4,cp4,'kv'); hold on; grid on;
        plot(lambdav,cp,'k-');  
        legend('single streamtube theory','multiple streamtube theory','Location','best');
        h=xlabel('$\lambda$','Interpreter','latex'); set(h,'Fontsize',20);
        h=ylabel('$C_P$','Interpreter','latex');
        set(h,'Fontsize',20);
        title('\sigma = 0.2 - Cd = 0');
end
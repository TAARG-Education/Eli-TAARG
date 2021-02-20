% main.m
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
% |Author      : Nasti Giuseppe, Tomasso Armando Diego                                         |
% |              University of Naples Federico II.                                             |
% |Version     : 1.1                                                                           |
% |Date        : 20/01/2021                                                                    |
% |Modified    : 18/02/2021                                                                    |
% |Description : test case for characteristic curves evaluation for a windmill                 |
%                                                                                              |
% |Reference   : Tognaccini R., (2019), "Lezioni di Aerodinamica dell'ala rotante"             |
%                BURTON, T., JENKINS, N.,SHARPE, D., BOSSANYI, E.,(2011),                      |
%                             "Wind Energy Handbook", 2nd Edition, Wiley and Sons.             |
% |Input       :  (c) = chord of the blade                                                     |
%                 (r) = radius station                                                         |
%                 (beta) = twist angle distribution                                            |
%                 (N) = number of blades                                                       |
%                 (First_lambda) = first value of tip speed                                    |
%                 (Last_lambda) = last value of tip speed                                      |
%                 (Aeromatrix)= matrix of aerodynamic characteristics                          |                                                            |
%                                                                                              |
% |Output      :  (Cp) = vector of power coefficient values for the                            |
%                 windmill                                                                     |
%                 (Cq) = vector of torque coefficient values for the                           |
%                 windmill                                                                     |
%                 (Ct) = vector of thrust coefficient values for the                           |
%                 windmill                                                                     |
%                 (Lambda) = vector of tip speed values corresponding to                       |
%                 the cp and cq coefficients                                                   |
% |Note        :                                                                               |
% ==============================================================================================

clc; clear all; close all;

N=3;                              %Number of Blades
test=load('BladeTestCase.txt');
R=20.5;                           %Radius of the blade [m]
r=test(:,1)*R;                    %radius stations [m]
c=test(:,2)*R;                    %Chord stations [m]
beta=test(:,3);                   %Beta stations [deg]
lambda_first=2;                   %First Lambda 
lambda_end=15;                    %Last Lambda

%Construction of the aerodynamic characteristics of blade stations 
Aero_matrix=ones(length(r),12);
Aero_matrix(:,1) = Aero_matrix(:,1)*0.0090;     %Cd_min
Aero_matrix(:,2) = Aero_matrix(:,2)*0.004;      %dCd_dCl2
Aero_matrix(:,3) = Aero_matrix(:,3)*0.3;        %Cl_Cd_min
Aero_matrix(:,4) = Aero_matrix(:,4)*4000000;    %Re_ref
Aero_matrix(:,5) = Aero_matrix(:,5)*4000000;    %Re_inf
Aero_matrix(:,6) = Aero_matrix(:,6)*(-0.2);     %f
Aero_matrix(:,7) = Aero_matrix(:,7)*1.5;        %Cl_max
Aero_matrix(:,8) = Aero_matrix(:,8)*(-0.78);    %Cl_min
Aero_matrix(:,9) = Aero_matrix(:,9)*(-3);       %alpha_zero_lift [deg]
Aero_matrix(:,10)= Aero_matrix(:,10)*0.11;      %Cl_alpha[1/deg]
Aero_matrix(:,11)= Aero_matrix(:,11)*0;         %Cl_alpha_stall [1/deg]

% Function calling
[Cp Cq Ct Lambda]=Characteristics_Curve_HO_Windmill(N,r,Aero_matrix,beta,c,lambda_first,lambda_end);
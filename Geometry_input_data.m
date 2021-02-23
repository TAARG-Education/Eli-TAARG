
%% \file Geometry Input
%  \brief starting from a database the function promyProp.v_infides different myProp.v_infariables about ...
...the propeller geometry that can be needed for different analysis 
%  \author Carmen Di Mauro & Teresa Aliperta
%  \myProp.v_infersion 1.0.0
%
% Eli-TAARG is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public
% License as published by the Free Software Foundation; either
% myProp.v_infersion 3 of the License, or (at your option) any later myProp.v_infersion.
%
% Eli-TAARG is demyProp.v_infeloped by the TAARG Educational organization for
% educational purposes only.
% Theoretical and Applied Aerodynamic Research Group - UnimyProp.v_infersity of Naples Federico II.
%
% Eli-TAARG GitHub link: <https://github.com/TAARG-Education/Eli-TAARG>
%
% Eli-TAARG is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without emyProp.v_infen the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
% General Public License for more details.
% <http://www.gnu.org/licenses/>.
%
% ==============================================================================================
% |Name        : Function for the propeller analysys input data                                 |
% |Author      : Carmen Di Mauro & Teresa Aliperta                                              |
% |              UnimyProp.v_infersity of Naples Federico II.                                              |
% |myProp.v_infersion     : 1.0.0                                                                          |
% |Date        : 25/11/2020                                                                     |
% |Modified    : 15/12/2020                                                                     |
% |Description : starting from a database the function promyProp.v_infides different myProp.v_infariables about the   | 
%                propeller geometry that can be needed for different analysis                   |
% |Reference   :                                                                                |
% |Input       : Function for the reading of the propeller geometry that from a file .txt reads |
%                the geometry of the propeller                                                  |
% |Output      : myProp.v_infariables that contain the geometry information needed                         |
% |Note        :                                                                                |
% ==============================================================================================

function [X]= Geometry_input_data(txt)

PropDataFileName='propgeometry.txt';    %.txt standard geometry file
myProp=geometryreader(PropDataFileName);       %class call


function_name=txt;
switch function_name
    case 'BEMT'
        X=struct('blades',myProp.N, 'r',myProp.r, 'hub_radius', myProp.h_r, 'velocity', myProp.v_inf, 'rpm',myProp.rpm);
    case 'Darrieus_flusso_multiplo'
        X=struct('corda', myProp.c, 'R', myProp.r, 'blades',myProp.N,'Cd', myProp.Cd);
    case 'Axial_rotor'
        X=struct('r',myProp.r);
    case 'Ang_attacco_effettimyProp.v_info'
        X=struct('theta',myProp.theta);
    case 'elica_intubata'
        X=struct('velocity',myProp.v_inf,'Thrust',myProp.T);
    case 'adim_coeff'
        X=struct('velocity',myProp.v_inf,'altitude',myProp.h,'blades', myProp.N);
    case 'Opti_Prop'
        X=struct('blades',myProp.N,'hub_radius',myProp.h_r,'tip_radius',myProp.t_r,'rpm',myProp.rpm,'altitude',myProp.h,'velocity', myProp.v_inf,'Thrust',myProp.T, 'Power', myProp.P);
    case 'RmyProp.v_infortexInt'
        X=struct('corda',myProp.c);
    case 'flappingangles'
        X=struct('velocity',myProp.v_inf, 'Lock number', myProp.lock_number,'R',myProp.r);
    case 'Cdcl_xfoil'
        X=struct('Reynolds number', myProp.Reynolds);
    case 'RotorFF'
        X=[]; %they only need in input the angle of attack
    case 'Axial_Descent_Ascent'
        X=struct('R',myProp.r);
    case 'Cdcl_xrotor'
        X=struct('Reynolds_number', myProp.Reynolds, 'f',myProp.f);
        
end



end



%% \file Geometry Input
%  \brief starting from a database the function provides different variables about ...
...the propeller geometry that can be needed for different analysis 
%  \author Carmen Di Mauro & Teresa Aliperta
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
% |Name        : Function for the propeller analysys input data                                 |
% |Author      : Carmen Di Mauro & Teresa Aliperta                                              |
% |              University of Naples Federico II.                                              |
% |Version     : 1.0.0                                                                          |
% |Date        : 25/11/2020                                                                     |
% |Modified    : 15/12/2020                                                                     |
% |Description : starting from a database the function provides different variables about the   | 
%                propeller geometry that can be needed for different analysis                   |
% |Reference   :                                                                                |
% |Input       : Function for the reading of the propeller geometry that from a file .txt reads |
%                the geometry of the propeller                                                  |
% |Output      : variables that contain the geometry information needed                         |
% |Note        :                                                                                |
% ==============================================================================================

function [X]= input_per_la_geometria(txt)
% PropDataFileName='propgeometry.txt';    %.txt standard geometry file
% myProp=geometryreader(PropDataFileName);       %class call
%      numdipale = myProp.N ;
    numdipale=3;
    diametro = 2;
    raggio=5;
    corda=5;
    pitch=[2 3 5 7];
    Reynolds=[ 10e5 20e6];
    Mach=5;
    rpm=3;
    hub_radius=5;
    tip_radius=5;
    v=50;
    n_sec=3;
    massa=500;
    theta=0;
    thrust=10;
    h=200;
    
function_name=txt;
switch function_name    
    case 'BEMT'
        X=struct('blades',numdipale, 'r',raggio, 'sections', n_sec, 'hub_radius', hub_radius, 'velocity', v, 'rpm',rpm);
    case 'Darrieus_flusso_multiplo'
        X=struct('pitch',pitch, 'Reynolds',Reynolds);
    case 'Axial_rotor'
        X=struct('mass',massa,'r',raggio);
    case 'Ang_attacco_effettivo'
        X=struct('theta',theta);
    case 'elica_intubata'
        X=struct('v',v,'thrust',thrust);
    case 'adim_coeff'
        X=struct('velocity',v,'altitude',h,'blades', numdipale);
    case 'Opti_Prop'
        X=struct('blades',numdipale,'hub_radius', hub_radius,'tip_radius',tip_radius,'rpm',rpm,'altitude',h,'velocity', v,'thrust',thrust);
        

end



end


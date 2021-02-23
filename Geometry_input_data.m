
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

function [X]= Geometry_input_data(txt)

PropDataFileName='propgeometry.txt';    %.txt standard geometry file
myProp=geometryreader(PropDataFileName);       %class call

numdipale = myProp.N ;
raggio= myProp.r ; %adimensional
corda=myProp.c; %adimensional
chords=myProp.chords; %distribution along the radius
radius=myProp.radius; %distribution along the radius
tip_radius= myProp.t_r ;
hub_radius=myProp.h_r ;
thrust=myProp.T;
power=myProp.P;
locknum=myProp.lock_number;
v=myProp.v_inf;
h=myProp.h;
theta=myProp.theta;
pitch=myProp.pitch;
Reynolds=myProp.Reynolds;
Mach=myProp.Mach;
Cd=myProp.Cd;
rpm=myProp.rpm;
f=myProp.f;
advanceratio=myProp.advanceratio;

function_name=txt;
switch function_name
    case 'BEMT'
        X=struct('blades',numdipale, 'r',raggio, 'hub_radius', hub_radius, 'velocity', v, 'rpm',rpm);
    case 'Darrieus_flusso_multiplo'
        X=struct('corda', corda, 'R', r, 'blades',numdipale,'Cd', Cd);
    case 'Axial_rotor'
        X=struct('r',raggio);
    case 'Ang_attacco_effettivo'
        X=struct('theta',theta);
    case 'elica_intubata'
        X=struct('v',v,'thrust',thrust);
    case 'adim_coeff'
        X=struct('velocity',v,'altitude',h,'blades', numdipale);
    case 'Opti_Prop'
        X=struct('blades',numdipale,'hub_radius', hub_radius,'tip_radius',tip_radius,'rpm',rpm,'altitude',h,'velocity', v,'thrust',thrust, 'power', power);
    case 'RvortexInt'
        X=struct('corda',corda);
    case 'flappingangles'
        X=struct('velocity',v, 'Lock number', locknum,'R',r);
    case 'Cdcl_xfoil'
        X=struct('Reynolds number', Reynolds);
    case 'RotorFF'
        X=[]; %they only need in input the angle of attack
    case 'Axial_Descent_Ascent'
        X=struct('R',r);
    case 'Cdcl_xrotor'
        X=struct('Reynolds_number', Reynolds, 'f',f);
        
end



end


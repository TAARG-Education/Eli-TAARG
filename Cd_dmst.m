%% \Cd_dmst.m
%  \brief: Returns Cd by interpolation or using a constant value of 100 DC
%  \author: Gabriele Lucci
%  \version: 1.1
%
% Eli-TAARG is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public
% License as published by the Free Software Foundation; either
% version 3 of the License, or (at your option) any later version.
%
% Eli-TAARG is developed by the TAARG Educational organization for
% educational purposes only.
% Theoretical and Applied Aerodynamic Research Group - University of Naples 
% Federico II.
%
% Eli-TAARG GitHub link: <https://github.com/TAARG-Education/Eli-TAARG>
%
% Eli-TAARG is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
% General Public License for more details.
% <http://www.gnu.org/licenses/>.
%
% =========================================================================
% |Name        : Cd.m                                                    
% |Author      : Gabriele Lucci                                            
% |              University of Naples Federico II.                         
% |Version     : 1.1                                                      
% |Date        : 02/09/2022                                                
% |Modified    : 02/16/2022                                                         
% |Description : Returns Cd value in two ways:
% |              1) through 2-variable interpolation over the data read and 
% |              provided by "ReadAeroData.m" function in the "RE", 
% |              "ALPHA", "Cd_data" global variables;
% |              2) using "ClCd_Xrotor.m" function in Eli-TAARG library;
% |Reference   : Sheldahl R. E. and Klimas P. C. (1981). "Aerodynamic 
% |              Characteristics of Seven Symmetrical Airfoil Sections 
% |              Through 180-Degree Angle of Attack for Use in Aerodynamic 
% |              Analysis of Vertical Axis Wind Turbines". Sandia National
% |              Laboratories.
% |Input       : (Re)       = double, Reynolds number;
% |              (alpha)    = double, angle of attack;
% |              (aeroflag) = string, 'xrotor' | 'skdata'. Choose between 
% |                           Cd calculation through "ClCd_Xrotor.m" 
% |                           ('xrotor') function included in Eli-TAARG 
% |                           library or through two-variable interpolation
% |                           on Sheldahl & Klimas (see references) 
% |                           experimental data ('skdata').
% |Output      : (Cd_val)   = drag coefficient value.
% |Note        : -
% =========================================================================
%
function Cd_val = Cd_dmst(Re,alpha,aeroflag)

global RE ALPHA Cd_data

if strcmpi(aeroflag,'skdata')
    
    Cd_val = interp2(RE,ALPHA,Cd_data,Re,alpha);
    
elseif strcmpi(aeroflag,'xrotor')    
    
    Cl = Cl_dmst(Re,alpha,aeroflag);
    input_v(1) = 0.005;
    input_v(2) = 0.0040;
    input_v(3) = 0;
    input_v(4) = 5e6;
    input_v(5) = Re;
    input_v(6) = -0.2;
    input_v(7) = 1.75;
    input_v(8) = -1.75;
    [~,~,Cd_val] = ClCd_XRotor(input_v, Cl);

else
    
    error("Spellcheck 'aeroflag'");
    
end

end
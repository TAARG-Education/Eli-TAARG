%% \Cd_dmst.m
%  \brief: Returns Cd by interpolation or using a constant value of 100 DC
%  \author: Gabriele Lucci
%  \version: 1.00
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
% |Version     : 1.00                                                      
% |Date        : 02/09/2022                                                
% |Modified    : -                                                         
% |Description : Returns Cd value in two ways:
% |              1) through 2-variable interpolation over the data read and 
% |              provided by "ReadAeroData.m" function in the "RE", 
% |              "ALPHA", "Cd_data" global variables;
% |              2) providing a constant value of 100 DC;
% |Reference   : Sheldahl R. E. and Klimas P. C. (1981). "Aerodynamic 
% |              Characteristics of Seven Symmetrical Airfoil Sections 
% |              Through 180-Degree Angle of Attack for Use in Aerodynamic 
% |              Analysis of Vertical Axis Wind Turbines". Sandia National
% |              Laboratories.
% |Input       : (Re)       = double, Reynolds number;
% |              (alpha)    = double, angle of attack;
% |              (aeroflag) = string, 'simple' | 'real'. Choose between 
% |                           Cl calculation through linear law (simple) 
% |                           or through two-variable interpolation
% |                           (realistic).
% |Output      : (Cd_val)   = drag coefficient value.
% |Note        : -
% =========================================================================
%
function Cd_val = Cd_dmst(Re,alpha,aeroflag)

global RE ALPHA Cd_data

if strcmpi(aeroflag,'real')
    
    Cd_val = interp2(RE,ALPHA,Cd_data,Re,alpha);
    
elseif strcmpi(aeroflag,'simple')

    Cd_val = 0.01;
    
else
    
    error("Spellcheck 'aeroflag'");
    
end

end
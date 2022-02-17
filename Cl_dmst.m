%% \Cl_dmst.m
%  \brief: Returns Cl by interpolation or using a linear law in alpha
%  \author: Gabriele Lucci
%  \version: 1.02
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
% |Name        : Cl.m                                                    
% |Author      : Gabriele Lucci                                            
% |              University of Naples Federico II.                         
% |Version     : 1.02                                                      
% |Date        : 02/09/2022                                                
% |Modified    : 02/17/2022                                                       
% |Description : Returns Cl value in two ways:
% |              1) through 2-variable interpolation over the data read and 
% |              provided by "ReadAeroData.m" function in the "RE", 
% |              "ALPHA", "Cl_data" global variables;
% |              2) employing a linear law in "alpha";
% |Reference   : Sheldahl R. E. and Klimas P. C. (1981). "Aerodynamic 
% |              Characteristics of Seven Symmetrical Airfoil Sections 
% |              Through 180-Degree Angle of Attack for Use in Aerodynamic 
% |              Analysis of Vertical Axis Wind Turbines". Sandia National
% |              Laboratories.
% |Input       : (Re)       = double, Reynolds number;
% |              (alpha)    = double, angle of attack;
% |              (aeroflag) = string, 'xrotor' | 'skdata' | 'simple'. 
% |                           Choose between Cl calculation through linear 
% |                           law ('xrotor') or through two-variable 
% |                           interpolation on Sheldahl & Klimas (see 
% |                           references) experimental data ('skdata' or
%                             'simple').
% |Output      : (Cl_val)   = lift coefficient value.
% |Note        : -
% =========================================================================
%
function Cl_val = Cl_dmst(Re,alpha,aeroflag)

global RE ALPHA Cl_data

if strcmpi(aeroflag,'skdata')
    
    Cl_val = interp2(RE,ALPHA,Cl_data,Re,alpha);
    
elseif strcmpi(aeroflag,'xrotor') || strcmpi(aeroflag,'simple')

    Cl_val = 2*pi*alpha;

else
    
    error("Spellcheck 'aeroflag'");
    
end

end
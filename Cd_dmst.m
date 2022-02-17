%% \Cd_dmst.m
%  \brief: Returns Cd by interpolation, with Xrotor eqn. or as a constant
%  \author: Gabriele Lucci
%  \version: 1.2
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
% |Version     : 1.2                                                      
% |Date        : 02/09/2022                                                
% |Modified    : 02/17/2022                                                         
% |Description : Returns Cd value in three ways:
% |              1) through 2-variable interpolation over the data read and 
% |              provided by "ReadAeroData.m" function in the "RE", 
% |              "ALPHA", "Cd_data" global variables;
% |              2) using "ClCd_Xrotor.m" function in Eli-TAARG library;
% |              3) as a constant equal to 100 DC.
% |Reference   : Sheldahl R. E. and Klimas P. C. (1981). "Aerodynamic 
% |              Characteristics of Seven Symmetrical Airfoil Sections 
% |              Through 180-Degree Angle of Attack for Use in Aerodynamic 
% |              Analysis of Vertical Axis Wind Turbines". Sandia National
% |              Laboratories.
% |Input       : (Re)       = double, Reynolds number;
% |              (alpha)    = double, angle of attack;
% |              (aeroflag) = string, 'xrotor' | 'skdata' | 'simple'. 
% |                           Choose between Cd calculation through 
% |                           "ClCd_Xrotor.m" ('xrotor') function included 
% |                           in Eli-TAARG library, through two-variable 
% |                           interpolation on Sheldahl & Klimas (see 
% |                           references) experimental data ('skdata'), or
% |                           setting it equal to 100 DC ('simple').
% |              (varargin) = 7-by-1 double array, input parameters for
% |                           "ClCd_XRotor.m" function (see documentation).
% |Output      : (Cd_val)   = drag coefficient value.
% |Note        : -
% =========================================================================
%
function Cd_val = Cd_dmst(Re,alpha,aeroflag,varargin)

narginchk(3,4);

if strcmpi(aeroflag,'xrotor') && nargin == 3

    error(['<strong>ClCd_XRotor</strong> also needs an input vector. ', ...
        'See function documentation.'])

end

global RE ALPHA Cd_data

if strcmpi(aeroflag,'skdata')
    
    Cd_val = interp2(RE,ALPHA,Cd_data,Re,alpha);

elseif strcmpi(aeroflag,'simple')
    
    Cd_val = 0.01;
    
elseif strcmpi(aeroflag,'xrotor')    
    
    Cl           = Cl_dmst(Re,alpha,aeroflag);
    input_v      = varargin{1};
    [~,~,Cd_val] = ClCd_XRotor(input_v,Re,Cl);

else
    
    error("Spellcheck 'aeroflag'");
    
end

end
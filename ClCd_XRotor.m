function [Cl_vec, Cd, Cd_bp] = ClCd_XRotor(input_v, Cl_bp)
%% \ClCd_XRotor.m
%  \brief: the function evaluates the Cd(Cl) law according to XRotor
%  aerodynamic model
%  \author: Crescenzo Visone
%  \version: 1.03
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
% |Name        : ClCd_XRotor.m                                                           |
% |Author      : Crescenzo Visone                                                        |
% |              University of Naples Federico II.                                       |
% |Version     : 1.03                                                                    |
% |Date        : 25/11/2020                                                              |
% |Modified    : 16/12/2020                                                              |
% |Description : Cd(Cl) evaluation law according to XRotor aero model                    |     
% |Reference   : XRotor User Guide                                                       |
% |              <http://web.mit.edu/drela/Public/web/xrotor/xrotor_doc.txt>             |                                                         |
% |Input       : * the first input must be a vector with the following components:       |
% |              (Cd_min) = minimum blade element's drag coefficient                     |
% |              (dCd_dCl2) = quadratic coefficient of the Cd(Cl) parable                |  
% |              (Cl_Cd_min) = lift coefficient @ minimum Cd                             |
% |              (Re_ref) = blade element's Reynolds number                              |
% |              (Re_inf) = asymptotic Reynolds number                                   |
% |              (f) = Reynolds number scaling exponent                                  |
% |              (Cl_max) = maximum blade element's lift coefficient                     |
% |              (Cl_min) = minimum blade element's lift coefficient                     |
% |              * the second input can either be a scalar value or a vector equal to :  | 
% |              (Cl_bp) = blade element's lift coefficient breakpoints                  |
% |Output      : (Cl_vec) = lift coefficient vector used to plot the law                 |
% |              (Cd) = drag coefficent values @ Cl_vec                                  |
% |              (Cd_bp) = drag coefficient values @ requested Cl_bp                     |
% |Note        :                                                                         |
% ==============================================================================================

    % taking input data
    Cd_min    = input_v(1);
    dCd_dCl2  = input_v(2);
    Cl_Cd_min = input_v(3);
    Re_ref    = input_v(4);
    Re_inf    = input_v(5);
    f         = input_v(6);
    Cl_max    = input_v(7);
    Cl_min    = input_v(8);
    % lift coefficient vector creation
    Cl_vec = Cl_min:.01:Cl_max;
    % aerodynamic model building
    Cd = (Cd_min + dCd_dCl2*(Cl_Cd_min - Cl_vec).^2)*(Re_inf/Re_ref)^f;
    %% evaluation of Cd values @ requested Cl breakpoints 
    % note that Cl_bp can be both a single value or a vector
    % definition of an anonymous function through the handle (@) symbol
    Cd_interp = @(Cl_interp) interp1(Cl_vec, Cd, Cl_interp, 'pchip');
    % anonymous fcn used to compute requested Cd
    Cd_bp = Cd_interp(Cl_bp);
    %% plots section ------------------------------------------------------
    % blade element's polar diagram
    % along x: drag coefficient
    % along y: lift coefficient
    plot(Cd, Cl_vec);
end
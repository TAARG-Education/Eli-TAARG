%% \main.m
%  \brief: this script tests the ClCd_XRotor.m function
%  \author: Crescenzo Visone
%  \version: 1.01
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
% |Name        : main.m                                                                  |
% |Author      : Crescenzo Visone                                                        |
% |              University of Naples Federico II.                                       |
% |Version     : 1.01                                                                    |
% |Date        : 22/12/2020                                                              |
% |Modified    : 22/12/2020                                                              |
% |Description : ClCd_XRotor.m function testing                    |     
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
% |Note        : results have been validated through XRotor software                     |
% ==============================================================================================
close all; clear all; clc;

% the following script is only aimed to test the ClCd_XRotor function

% input vector
input_vec = [.0068, .0023, .69, 750000, 750000, -1.5, 1.57, -.86];

% lift coefficient input breakpoints
Cl_bp = [0.8,1];

% function calling
[Cl_vec, Cd, Cd_bp] = ClCd_XRotor(input_vec, Cl_bp);
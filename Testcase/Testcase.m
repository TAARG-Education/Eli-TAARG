%% \main.m
%  \brief  Graphs to describe propellers, rotors and turbines operation
%  \author Giada Cardellino & Antonio Golino
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
% |Name        :                                                                               |
% |Author      : Giada Cardellino & Antonio Golino                                             |
% |              University of Naples Federico II.                                             |
% |Version     : 1.0.0                                                                         |
% |Date        : 25/11/20                                                                      |
% |Modified    : 27/12/20                                                          	           |
% |Description : starting from the output data of all functions, the function		           |
% |              provides to show in output all the charateristics graphs to describe          |
% |              propellers, rotors and turbines operation 				                       |
% |											                                                   |
% |Reference   :                                                                               |
% |Input       : output data of all functions                                                  |
% |Output      : graphs to describes propellers, rotors and turbines operation     	           |
% |                                                                                            |
% |Note        :                                                               		           |
% ==============================================================================================
close all; clear all; clc;

% the following script is only created to test the Output Function

% input for the function: these vectors are randomly invented to test the function but there is no case ClCd_Xrotor reference 
cd=linspace(0,1000,1000);
cl=linspace(-1,1,1000);

% calling the output function to give a plot
[y1,y2]=FunzionidiOutput('ClCd_XRotor',cd,cl);
% calling the output function to give a text file
[y1,y2]=FunzionidiOutput('Opti_prop_T');



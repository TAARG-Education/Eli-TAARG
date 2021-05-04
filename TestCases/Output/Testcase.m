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
% |		 propellers, rotors and turbines operation 				                               |
% |											                                                   |
% |Reference   :                                                                               |
% |Input       : 'Function name' and vectors to plot in this way:                              |
%                  v1=x-axis, v2=y-axis (plot 1)                                               |
%                  v3=x-axis, v4=y-axis (plot 2)                                               |
%                  (if there isn't a plot 2 v3=[]; v4=[]                                       |
%                 'Function name if the output is a data text file with                        |
%                 vi=[]                                                                        |
% |Output      : graphs to describes propellers, rotors and turbines operation or              |
%                data text file                                                                |
% |                                                                                            |
% |Note        :                                                               		           |
% ==============================================================================================
close all; clc;

% the following script is only created to test the Output Function

% input for the function: these vectors are randomly invented to test the function but there is no case ClCd_Xrotor reference 
cl=matrix(:,2);
cd=matrix(:,3);
% calling the output function to give a plot
[y1,y2,y3]=Outputfunction('ClCd_XRotor',cd,cl,[],[],[],[]);




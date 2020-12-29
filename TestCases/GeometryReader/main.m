%% \file funzione per la lettura della geometria
%  \brief the function reads the geometry of a propeller and analyzes it   
%         using xrotor
%  \author Davide Mennella & Alessandro Carmine Esposito
%  \version 1.0
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
% |Name        : Function for the reading of the propeller geometry                            |
% |Author      : Davide Mennella& Alessandro Carmine Esposito                                  |
% |              University of Naples Federico II.                                             |
% |Version     : 1.0                                                                           |
% |Date        : 16/11/2020                                                                    |
% |Modified    : 29/12/2020                                                                    |
% |Description : the function reads the geometry of a propeller and analyzes it using xrotor   |
% |Reference   : http://web.mit.edu/drela/Public/web/xrotor/                                   |
% |Input       : standard .txt file that contains informations about propeller geometry        |
% |Output      : local variables that contain the geometry information about the propeller     |
%                and/or propeller analysis via interface with xrotor software                  | 
% |Note        : for the correct use of the tool, user has to put in the same folder:          |
%                - the matlab file 'main';                                                     |
%                - the class file 'geometryreader';                                            |
%                - the xrotor executable;                                                      |
%                - the standard 'propgeometry.txt' file with the geometry parameters           |           
% ==============================================================================================


%% User instructions analize and plot global variables must be setted as follow:

%  The Matlab class 'geometryreader' reads parameters from the input file 
%  'propgeometry.txt' and allocates them in local variables.
%  Moreover the class itself contains two functions:
%   - 'xRotorAnalysis' that gives the propeller's variables in input to 
%     xrotor.exe which analizes the propeller in the range of advance ratio 
%     setted in the propgeometry.txt and saves the analysis  results in 
%     local variables. User has to call this function with the following
%     sintax -> myProp.xRotorAnalysis after calling the class.
%   - 'plotResults' that performs the analysis and also plots the analysis
%     results. User has to call this function with the following sintax ->
%     myProp.plotResults after calling the class.

%% main
clc; clear; close all;
PropDataFileName='propgeometry.txt';    %.txt standard geometry file
myProp=geometryreader(PropDataFileName);       %class call
%myProp.xRotorAnalysis;
%myProp.plotResults;



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
% |Modified    : 20/11/2020                                                                    |
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



%The analize and plot global variables must be setted as follow:

%- analize=0: the class 'geometryreader' only reads parameters from the 
%  input file 'PropDataFileName' and allocates them in local variables

%- analize=1 && plots=0: the class 'geometryreader' reads parameters from 
%  the input file 'PropDataFileName' and allocates them in local variables; 
%  moreover the function give the propeller's variables in unput to 
%  xrotor.exe which analizes the propeller in the range of advance ratio 
%  setted in the variables firstadv, lastadv,advstep and saves the 
%  analysis  results in local variables;

%- analize=1 && plots=1: the class also generates plots of the analysis 
%  results.

clc; clear; close all;
global analize plots firstadv lastadv advstep
firstadv='0.143';  %first value of the advance ratio vector of the analysis
lastadv='0.363';   %last value of the advance ratio vector of the analysis
advstep='0.003';   %step of the advance ratio
analize=1;
plots=0;
PropDataFileName='propgeometry.txt';    %.txt standard geometry file
myProp=analisi(PropDataFileName);       %class call


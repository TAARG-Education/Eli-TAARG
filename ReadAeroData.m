%% \mst.m
%  \brief: Read NACA 0012 aerodynamic force coefficients from file
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
% |Name        : ReadAeroData.m                                                    
% |Author      : Gabriele Lucci                                            
% |              University of Naples Federico II.                         
% |Version     : 1.00                                                      
% |Date        : 02/09/2022                                                
% |Modified    : -                                                         
% |Description : Reads aerodynamic data stored in a Microsoft Excel
% |              spreadsheet. The data is stored in global variables to be 
% |              recalled by lift and drag interpolation routines, 
% |              toghether with angle of attack and Reynolds number grids.
% |Reference   : Sheldahl R. E. and Klimas P. C. (1981). "Aerodynamic 
% |              Characteristics of Seven Symmetrical Airfoil Sections 
% |              Through 180-Degree Angle of Attack for Use in Aerodynamic 
% |              Analysis of Vertical Axis Wind Turbines". Sandia National
% |              Laboratories.
% |Input       : (filepath) = string, data spreadsheet path.
% |Output      : -
% |Note        : Aerodynamic data spreadsheet must be formatted as in
% |              "sandia0012data.xlsx".
% =========================================================================
%
%%
function ReadAeroData(filepath)

global Cl_data Cd_data ALPHA RE RADrunflag

RADrunflag = true;

aerodata = readtable(filepath);

Re_data = aerodata{~isnan(aerodata{:,end}),end};

for i = 1:11
    
    j = 2*i;
    Cl_data(:,i) = aerodata{:,j};
    
end

for i = 1:11
    
    j = (2*i+1);
    Cd_data(:,i) = aerodata{:,j};
    
end

alpha = [flipud(-aerodata.alpha(2:end)); aerodata.alpha];

Cl_data = [flipud(-Cl_data(2:end,:));Cl_data];
Cd_data = [flipud(Cd_data(2:end,:));Cd_data];

[RE,ALPHA] = meshgrid(Re_data,alpha);

end
%% \xfoil_polar.m
%  \brief: evaluation of the Cd(Cl) polar through the software xfoil
%  \author(s): Palma Caputo - Lorenzo Frascino
%  \version: 1.5
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
% |Name        : CdCl_xfoil.m                                                                    |
% |Author(s)   : Palma Caputo - Lorenzo Frascino                                                 |
% |              University of Naples Federico II.                                               |
% |Version     : 1.5                                                                             |
% |Date        : 03/01/2020                                                                      |
% |Modified    : 03/01/2020                                                                      |
% |Description : evaluation of the Cd(Cl) polar by using the software xfoil                      |     
% |Reference   : xfoil user's Guide                                                              |
% |              <https://web.mit.edu/drela/Public/web/xfoil/xfoil_doc.txt>                      |                                                         
% |Input       : (airfoil) = it is possible to insert a NACA number (ie. 2412)                   |
% |                          or an external geometry file (ie. 'myairfoil.txt')                  |
% |              (numPanel) = number of panels for the geometry   (ie. 160)                      |
% |              (Re_number) = Reynolds number (ie. 1e6)                                         |
% |              (FirstAlfa) = first value of the angle of attack (deg) (ie. -10)                |
% |              (LastAlfa) = last value of the angle of attack (deg)   (ie. 10)                 |
% |              (DeltaAlfa) = pace between the angles of attack        (ie. 1)                  |
% |Output      : (Cl) = lift coefficient vector                                                  |
% |              (Cd) = drag coefficent vector                                                   |
% |Note        : The function must be in the same directory of xfoil.exe!                        |
% ==============================================================================================
function [Cl, Cd] = CdCl_xfoil(airfoil, numPanel, Re_number, FirstAlfa, LastAlfa, DeltaAlfa)

numPanel_st  = num2str(numPanel);
Re_number_st = num2str(Re_number);
FirstAlfa_st = num2str(FirstAlfa);
LastAlfa_st  = num2str(LastAlfa);
DeltaAlfa_st = num2str(DeltaAlfa);

iter = '100';
Alfa_vec = FirstAlfa:DeltaAlfa:LastAlfa;

saveGeometry = 'Airfoil_geometry.txt';  % Create .txt file to save airfoil coordinates
savePolar    = 'Polar.txt';             % Create .txt file to save the polar


% Delete files if they exist
if (exist(saveGeometry,'file'))
    delete(saveGeometry);
end

if (exist(savePolar,'file'))
    delete(savePolar);
end


%% WRITING XFOIL COMMANDS
% Create the airfoil
f_input = fopen('xfoil_input.txt','w');            % Create input file for xfoil 
fprintf(f_input,'y\n');

if isnumeric(airfoil)
    airfoil_1 = num2str(airfoil);
    fprintf(f_input,['naca ' airfoil_1 '\n']);
else
    fprintf(f_input,['load ' airfoil '\n']);
end

fprintf(f_input,'PPAR\n');                          
fprintf(f_input,['N ' numPanel_st '\n']);
fprintf(f_input,'\n\n');

% Data for the polar
fprintf(f_input,'OPER\n');                        
fprintf(f_input,'visc\n');                            
fprintf(f_input,[Re_number_st '\n']);                  
fprintf(f_input,['iter ' iter '\n']);
fprintf(f_input,'pacc\n');
fprintf(f_input,[savePolar '\n\n']);

for i = 1:length(Alfa_vec)
    fprintf(f_input,'a %2.4f\n', Alfa_vec(i));
end

fprintf(f_input,'pacc\n\n');

% Save the airfoil data points
fprintf(f_input,['PSAV ' saveGeometry '\n']);

% Close file
fclose(f_input);


%% RUNNING XFOIL (MUST BE IN THE SAME DIRECTORY!)
cmd = 'xfoil.exe < xfoil_input.txt';
[status,result] = system(cmd);

%% READ DATA FILE
filePol = fopen(savePolar);
A = textscan(filePol,'%f %f %f %f %f %f %f', 'Headerlines',12);
fclose(filePol);
alfa = A{1}(:,1);
Cl   = A{2}(:,1);
Cd   = A{3}(:,1);

figure(1);
plot(Cd,Cl,'k.-')
xlabel('Drag coefficient C_d');
ylabel('Lift coefficient C_l');
grid on;
end

%% \file nomefile
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
% |Description : starting from the output data of all functions, the function		               |
% |              provides to show in output all the charateristics graphs to describe          |
% |		 propellers, rotors and turbines operation 				                                       |
% |											                                                                       |
% |Reference   :                                                                               |
% |Input       : output data of all functions                                                  |
% |Output      : graphs to describes propellers, rotors and turbines operation     	           |
% |                                                                                            |
% |Note        :                                                               		             |
% ==============================================================================================

function [y1,y2] = Outputfunction(txt,v1,v2)

%% Default case
% function y = FunzionidiOutput(input1,input2,axisname)
% [rows,columns]=size(input1);
% for i=1:rows
%     xaxisname=string(axisname(1,i));
%     yaxisname=string(axisname(2,i));
%     y(i)=figure(i);
%     plot(x(i,:),y(i,:),'-k');
%     grid on;
%     xlabel(xaxisname);
%     ylabel(yaxisname);
% end
% end

%% Specifics cases
function_name=txt;
switch function_name
case 'ClCd_XRotor'
    y1=figure(1)
    plot(v1,v2,'linewidth',1.1);
    grid on;
    xlabel('C_d');
    ylabel('C_l');
    y2=[];
    case 'Axial_Descent_Ascent'
    y1=figure(1)
    plot(v1,v2,'linewidth',1.1);
    grid on;
    xlabel('$\widetilde{V}$','Interpreter','latex','FontSize',15);
    ylabel('$\widetilde{w}$','Interpreter','latex','FontSize',15);
    y2=figure(2)
    plot(v1,v2,'linewidth',1.1);
    grid on;
    xlabel('$\widetilde{V}$','Interpreter','latex','FontSize',15);
    ylabel('$\widetilde{P}$','Interpreter','latex','FontSize',15);
    case 'RVortexInt'
    y1=figure(1)
    plot(v1,v2);
    grid on;
    axis ([0 1 0 1]);
    text(0.25,1,'Velocity induced by vortex ring:');
    text(0.25,0.90,['fx= ',num2str(fx)]);
    y2=[];
    case 'CdCl_xfoil'
    y1=figure(1)
    plot(v1,v2);
    grid on;
    xlabel('Drag coefficient C_d');
    ylabel('Lift coefficient C_l');
    grid on;
    y2=[];
    case 'RotorFF' %subplot
    case 'Opti_prop_P'
    y1 = ['Data_Opti_Prop_P.txt'];
    fid = fopen(y1, 'wt');
    fprintf(fid, '%s\t%s','  efficiency =',eta,' at J =', J);  % header
    fprintf(fid,'\n');
    fprintf(fid, '%s\t%s','  w_conv =',w);  % header
    fprintf(fid,'\n');
    fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\n', '  r_adim','  chi','      a(chi)','  a''(chi)','  dCt/dr_adim',' dCp/dr_adim');  % header
    fclose(fid);
    dlmwrite(y1,DATA,'delimiter','\t','precision',['%10.',num2str(6),'f'],'-append');
    y2=[];
    case 'Opti_prop_T'
    y1 = ['Data_Opti_Prop_T.txt'];
    fid = fopen(y1, 'wt');
    fprintf(fid, '%s\t%s','  efficiency =',eta,' at J =', J);  % header
    fprintf(fid,'\n');
    fprintf(fid, '%s\t%s','  w_conv =',w);  % header
    fprintf(fid,'\n');
    fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\n', '  r_adim','  chi','      a(chi)','  a''(chi)','  dCt/dr_adim',' dCp/dr_adim');  % header
    fclose(fid);
    dlmwrite(y1,DATA,'delimiter','\t','precision',['%10.',num2str(6),'f'],'-append');  
    y2=[];
end
end




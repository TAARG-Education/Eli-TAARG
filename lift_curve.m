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
% |Name        : lift_curve.m                                                                  |
% |Author      : Nasti Giuseppe, Tomasso Armando Diego                                         |
% |              University of Naples Federico II.                                             |
% |Version     : 1.1                                                                           |
% |Date        : 20/01/2021                                                                    |
% |Modified    : 14/02/2021                                                                    |
% |Description : This function approximates lift curve with three segments:                    |
% |                - linear part of the lift-curve                                             |
% |                - post stall sector with Cl>0                                               |
% |                - post stall sector with Cl<0                                               |
% |                                                                                            |
% |Reference   :                                                                               |
% |Input       : the first input must be a vector with the following components:               |
% |              - cl_max                                                                      |
% |              - cl_min                                                                      |
% |              - alpha_zero_lift [deg]                                                       |
% |              - cl_alpha;       [1/deg]                                                     |
% |              - dcldalfastall   [1/deg] : value of Cl_alpha post stall                      |
% |                                                                                            |
% |              angle [deg]: Value of alpha where C_l is evaluate                             |
% |Output      : CL: Value of Cl(angle)                                                        |
% |Note        :                                                                               |
% ==============================================================================================


function CL=lift_curve(v_input,angle)

%Input mangament
cl_max=v_input(1);
cl_min=v_input(2);
alpha_zero_lift=v_input(3);     %deg
cl_alpha=v_input(4);            %1/deg
dcldalfastall=v_input(5);       %1/deg

v_alpha=linspace(-60,60,300);

%Ideal section of the lift curve
cl=cl_alpha*(v_alpha-alpha_zero_lift);

alphastallmax=(cl_max/cl_alpha)+alpha_zero_lift; % Ideal case
alphastallmin=(cl_min/cl_alpha)+alpha_zero_lift; %Ideal case

for i=1:length(cl)
    
    if cl(i)>cl_max %Correction for post stall with Cl>0

    cl(i)=cl_max+dcldalfastall*(v_alpha(i)-alphastallmax);
        
    elseif cl(i)<cl_min %Correction for post stall with Cl<0
        
    cl(i)=cl_min+dcldalfastall*(v_alpha(i)-alphastallmin);
        
    end
    
end

%Interpolation
CL=interp1(v_alpha,cl,angle,'linear');

end
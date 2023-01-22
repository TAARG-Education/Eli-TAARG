%% \RateOfClimb.m
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
% |Name        : RateOfClimb                                                                    |
% |Author      : Rosa Castiello & Raffaella Scarano                                             |
% |              University of Naples Federico II                                               |
% |Version     : 1.0.0                                                                          |
% |Date        : 07/01/2023                                                                     |  
% |Description : the function calculates rate of climb at different quotes and weights.         |        
% |                                                                                             |
% |Reference   : "Lezioni di Aerodinamica dell'ala rotante" by Renato Tognaccini.               |
% |Input       : Helicopter Weigth W, Total Required Power, Available Power,                    |
% |              Asymptotic Current Speed.                                                      |
% |                                                                                             |  
% |Output      : Rate Of Climb, Fastest Climb Speed, Maximum Rate Of Climb                      |
% |              Steepest Climb Speed, Angle of Maximum Rate Of Climb                           |
% |              Speed of maximum endurance and Speed of maximum Range.                         | 
% |Note        :                                                                                |                                               
% ==============================================================================================


function [RC,...
          Vy, RC_rapid,...
          Vx, gamma_ripid] = RateOfClimb(Wto, Ptot_req,Pavai,Vinf)


     RC = (Pavai - Ptot_req)./Wto;   %(Pavai-Ptot_req) is Power in excess
    
    for i=1:length(RC)
        if RC(i)<0         
            RC(i) = 0;   % If RC<0 we don't have the Required Power to climb
        end
    end
   

    for i=1:length(Vinf)
       if (RC(i) == max(RC))
           Vy = Vinf(i);    % Vy is the fastest climb speed, 
                            %at which we have the maximum rate of climb  
           RC_rapid = RC(i);
       end
       if (RC(i)/Vinf(i) == max(RC./Vinf))
           Vx = Vinf(i);              % Vx is the steepest climb, at which 
                                      %we have the maximum climbing angle
           RC_ripid = RC(i);
           gamma_ripid = asin(RC_ripid/Vx)*180/pi;
          
       end
    end

end
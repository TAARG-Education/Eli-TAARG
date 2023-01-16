%% \PowersInForwardFlight_main_tail_rotor.m
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
% |Name        : PowersInForwardFlight_main_tail_rotor                                          |
% |Author      : Rosa Castiello & Raffaella Scarano                                             |
% |              University of Naples Federico II                                               |
% |Version     : 1.0.0                                                                          |
% |Date        : 07/01/2023                                                                     |  
% |Description : the function calculates the required power in forward flight                   |
% |              for main rotor and tail rotor, according to Rotor BEMT theory,                 |        
% |              at different quotes and weights.                                               |
% |Reference   : "Lezioni di Aerodinamica dell'ala rotante" by Renato Tognaccini.               |
% |Input       : Altitude h, velocity vector, helicopter weigth W, Main rotor's radius,         |
% |              Main rotor's Area, Main rotor's Solidity, Main rotor's angular velocity,       |
% |              Cd of main rotor's airfoil, Tail rotor's radius, Tail rotor's Area,            |
% |              Tail rotor's solidity, Tail rotor's angular velocity,Cd of tail rotor's airfoil|
% |              b, correction factor, Auxiliary Power, Coefficient due to transmission losses, |
% |              Available Power at sea level.                                                  |  
% |Output      : Induced Power, Parasitic Power and Required Power by the main rotor, Parasitic |
% |              Power by the fuselage, Induced Power, Parasitic Power and Required Power       |
% |              by the tail rotor, Total Required Power, Available Power, Maximum Speed,       |
% |              Speed of maximum endurance and Speed of maximum Range.                         | 
% |Note        :                                                                                |                                               
% ==============================================================================================
%Argument 
% v_input:
% v_input(1) --> Main rotor's radius
% v_input(2) --> Main rotor's Area
% v_input(3) --> Main rotor's Solidity
% v_input(4) --> Main rotor's angular velocity
% v_input(5) --> Cd of main rotor's blade element
% v_input(6) --> Tail rotor's radius
% v_input(7) --> Tail rotor's Area
% v_input(8) --> Tail rotor's Solidity
% v_input(9) --> Tail rotor's angular velocity
% v_input(10)--> Cd of tail rotor's blade element
% v_input(11)--> b: distance between main rotor and tail rotor
% v_input(12)--> correction factor
% v_input(13)--> Auxiliary Power
% v_input(14)--> Coefficient due to transmission losses 
% v_input(15)--> Available Power at sea level

function [ Pind_B, Ppar_rp,Pfus,P_tot_rp,...
         Pind_rc,Ppar_rc,P_tot_rc,Ptot_req, Pdisp, ...
         V_NE, V_BE, V_BR] = PowerInForwardFlight_main_tail_rotor(h, Vinf, W,v_input)


% Required Power by the main rotor in forward flight
[~, ~, ~, rho]     = atmosisa(h);  
[~, ~, ~, rho_sl]  = atmosisa(0);

Pdisp = v_input(15)*rho/rho_sl;

%% Total power of the main rotor
%Parasitic Power
mu_rp   = Vinf./(v_input(4)*v_input(1));
Ppar_rp = (v_input(3)*v_input(5)/8.*(1+4.7*mu_rp.^2)).*...
          (pi*rho*v_input(4)^3*v_input(1)^5); % Rotor BEMT
                                      % obs: k=4.7 for radial speeds


T = W;  % Hp small alpha, in forward flight
%Induced Power 

%Impulsive Theory
Pind_A     = W^2./(2*rho*v_input(2).*Vinf); %Glauert Theory

w_hov      = 1/sqrt(2*rho)*sqrt(T/v_input(2));
Vinf_tilde = Vinf./w_hov;
w_tilde    = sqrt((- (Vinf_tilde.^2)./2 + sqrt(((Vinf_tilde.^2)./2).^2+1)));
w_rotore   = w_tilde*w_hov;
Pi_tilde   = w_tilde; %Pi_tilde = Vinf_tilde*sin(alpha) + w = w because sin(alpha)=0
                      % Hp: T costant
P_hov  = T*w_hov;

Pind_B = Pi_tilde*P_hov;

%Fusulage Power
f    = 0.007*v_input(2);         % f: equivalent wet area
Pfus = 0.5*rho.*Vinf.^3*f;  
%real main rotor's request Power
P_tot_rp = v_input(12)*Pind_B + Ppar_rp + Pfus;

%% Total Power of the tail rotor 

mu_rc = Vinf./(v_input(9)*v_input(6));

% Induced Power
Q_rot_main      = P_tot_rp./(v_input(4));
T_tail_required = Q_rot_main/v_input(11);
whov_rc         = 1/sqrt(2*rho).*sqrt(T_tail_required./v_input(7)); % Impulsive Theory
Vinf_tilde_rc   = Vinf./whov_rc;

w_tilde_rc  = sqrt((- (Vinf_tilde_rc.^2)./2 +...
              sqrt(((Vinf_tilde_rc.^2)./2).^2+1)));
w_rotor_rc  = w_tilde_rc.*whov_rc;

Pind_rc     = T_tail_required.*w_rotor_rc;

%Parasitic Power
Ppar_rc  = (v_input(8)*v_input(10)/8.*(1+4.7*mu_rc.^2)).*...
            pi*rho*v_input(9)^3*v_input(6)^5;
% Total power of the tail rotor
P_tot_rc = v_input(12)*Pind_rc + Ppar_rc;

% Total power required by total elicopter system
% in forward flight

v_P_aus  = v_input(13)*ones(1,length(Vinf));
Ptot_req = (P_tot_rp+P_tot_rc+v_P_aus)*v_input(14);

% Maximum speed in forward flight
for f=(length(Vinf)/2):length(Vinf)
  if(Pdisp>=Ptot_req(f-1) && Pdisp<Ptot_req(f))
       V_NE = (Vinf(f)+Vinf(f-1))/2;
  end
end

% Speed of maximum endurance
for f=1:length(Vinf)
 if (Ptot_req(f)==min(Ptot_req))
       V_BE = Vinf(f);
 end
end

% Speed of maximum Range
for f=1:length(Vinf)
 if(Ptot_req(f)/Vinf(f) == min(Ptot_req./Vinf))
      V_BR = Vinf(f);
 end
end

end


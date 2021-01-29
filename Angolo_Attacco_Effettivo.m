%% \file Angolo_Attacco_Effettivo
%  \brief Calculation and visualization of the blade section angle of
%  attack during a complete rotation of the rotor blade. 
% \author Serena Manganiello, Michele Massa. 
% \version 1.0
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
% |Name        : Angolo_Attacco_Effettivo
% |Author      : Serena Manganiello, Michele Massa.
% |              University of Naples Federico II.
% |Version     : 1.0                                                                      
% |Date        : 01/12/2020                                                            
% |Modified    : 16/01/2021                                                        
% |Description : Calculation and visualization of the effective angle of attack
%                during rotation of the rotor blade.                                             
% |Reference   : R. Tognaccini. "Lezioni di aerodinamica dell’ala rotante" 2019 pp. 109-112.
% |Input       : < lambda > rotor inflow ratio.
%                < r_segn > adimensional coordinate along the rotor blade.
%                < beta > blade flap angle.
%                < dbeta > blade flap angle's derivative.
%                < mu > rotor advance ratio.
%                < psi > azimuth angle of the blade.
%                < theta > blade pitch.
%                < alpha_stall_up > maximum blade section angle of attack.
%                < alpha_stall_lo > minimum blade section angle of attack.
% |Output      : blade section angle of attack.
% |Note        : u_R was calculated to be thorough; stall angles are
%                imposed in order to visualize only the non-stalled region
%                of the rotor.
% ==============================================================================================

function [alpha_e] = Angolo_Attacco_Effettivo(lambda,r_segn,beta,dbeta,mu,psi,theta,alpha_stall_up,alpha_stall_lo)

% Creation of the mesh grid.
[r_2d,psi_2d] = meshgrid(linspace(0,r_segn(1),length(r_segn)),psi);
x_hub = r_2d.*cos(psi_2d-pi/2);
y_hub = r_2d.*sin(psi_2d-pi/2);

[r_2d,psi_2d] = meshgrid(r_segn,psi);
x = r_2d.*cos(psi_2d-pi/2);
y = r_2d.*sin(psi_2d-pi/2);

% Variables' initialization.
u_P = zeros(length(r_segn),length(psi)); % air velocity of the blade section, perpendicular to the disk plane.
u_T = zeros(length(r_segn),length(psi)); % air velocity of the blade section, tangent to the disk plane.
u_R = zeros(length(r_segn),length(psi)); % radial air velocity of blade section.
phi = zeros(length(r_segn),length(psi)); % section inflow angle.
alpha_e = zeros(length(r_segn),length(psi)); % blade section angle of attack.

for i = 1:length(psi)
    u_P(:,i) = lambda + r_segn'.*dbeta(i) + beta(i).*mu.*cos(psi(i));
    u_T(:,i) = r_segn' + mu.*sin(psi(i));
    u_R(:,i) = mu*cos(psi(i));
    phi(:,i) = u_P(:,i)./u_T(:,i); 
    alpha_e(:,i) = theta' - phi(:,i);
end

stall = zeros(length(r_segn),length(psi)); % stalled region of the rotor blade.
non_stall = NaN(length(r_segn),length(psi)); % non-stalled region of the rotor blade.
for iii = 1:length(r_segn)
    for jjj = 1:length(psi)
        if alpha_e(iii,jjj) - convang(alpha_stall_up,'deg','rad') > 0 || alpha_e(iii,jjj) - convang(alpha_stall_lo,'deg','rad') < 0
            stall(iii,jjj) = alpha_e(iii,jjj);
        else
            non_stall(iii,jjj) = alpha_e(iii,jjj);
        end
    end
end

% Plot.
figure(1)
set(figure(1),'Color','w');
clf
fill(x_hub,y_hub,'k'); hold on; axis equal; axis off; grid off
text(r_segn(end-10),r_segn(end),['\mu = ',num2str(mu)],'Color','r','FontSize',12);
contourf(x,y,convang(non_stall','rad','deg'),20,'ShowText','on'); drawnow; axis equal;
end
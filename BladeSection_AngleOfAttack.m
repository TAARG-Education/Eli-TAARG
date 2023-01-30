%% \file Angolo_Attacco_Effettivo
% \brief Calculation and visualization of the blade section angle of
%        attack and effective Mach number during a complete rotation 
%        of the rotor blade.
% \author Serena Manganiello, Michele Massa.
% \updated Gianluca Porpora.
% \version 2.0
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
% |Author      : Serena Manganiello, Michele Massa, Gianluca Porpora.
% |Updated     : Gianluca Porpora
% |              University of Naples Federico II.
% |Version     : 2.0
% |Date        : 01/12/2020
% |Modified    : 30/01/2023
% |Description : Calculation and visualization of the effective angle of attack
%                and effective Mach during rotation of the rotor blade.
%                Function also studies 2D static stall and plots stalled
%                sections of the rotor blade. The solution of the
%                articulated rotor is required.
% |Reference   : R. Tognaccini. "Lezioni di aerodinamica dell’ala rotante" 2019 pp. 109-112.
% |Input       : < r_segn > adimensional coordinate along the rotor blade
%                           starting for the hub station.
%                < psi > azimuth angle of the blade (rad).
%                < mu > rotor advance ratio calculated as
%                       (V_inf*cos(alpha_tpp))/(Omega_tip)
%                < alpha_tpp > angle between V_inf and blade tip plane (deg)
%                < lambda > rotor inflow ratio.
%                < beta > blade flap angle (rad).
%                < dbeta > blade flap angle's derivative (rad).
%                < theta > blade pitch (rad)
%                < h > altitude for speed of sound (mt)
%                < tip_speed > for Mach number expressed in mt/sec
%                < alpha_stall_up > positive stall angle of attack.
%                < alpha_stall_lo > negative stall angle of attack.
%                < color > plot settings. Function accepts 'on' for contour
%                          plots and 'off' for black isolines.
% |Output      : blade section angle of attack, Mach number and graphics
% |Note        : u_R was calculated to be thorough; stall angles are
%                imposed in order to visualize only the non-stalled region
%                of the rotor. If there is at least one stalled section
%                function return a warning and a representation of the
%                stalled area of the disk.
% ==============================================================================================

function [alpha_eff,Meff] = BladeSection_AngleOfAttack(r_segn,psi,mu,alpha_tpp,lambda, ...
    beta,dbeta,theta,h,tip_speed,alpha_stall_up,alpha_stall_lo,color)

%% Calculate effective AoA, Mach and stall region

% Variables' initialization.
u_P = zeros(length(r_segn),length(psi));        % air velocity of the blade section, perpendicular to the disk plane.
u_T = zeros(length(r_segn),length(psi));        % air velocity of the blade section, tangent to the disk plane.
u_R = zeros(length(r_segn),length(psi));        % radial air velocity of blade section.
phi = zeros(length(r_segn),length(psi));        % section inflow angle.
alpha_e = zeros(length(r_segn),length(psi));    % blade section angle of attack.
Meff = zeros(length(r_segn),length(psi));       % Effective Mach
% Sound velocity calculation
[~,sound_speed,~,~] = atmosisa(h);

% Compute and distinguish stalled stations and inverse flow region
for i = 1:length(psi)
    u_P(:,i) = lambda + r_segn'.*dbeta(i) + beta(i).*mu.*cos(psi(i));
    u_T(:,i) = r_segn' + mu.*sin(psi(i));
    u_R(:,i) = mu*cos(psi(i));
    phi(:,i) = u_P(:,i)./u_T(:,i);
    alpha_e(:,i) = theta' - phi(:,i);
    Meff(:,i) = ((sqrt(u_P(:,i).^2 + u_T(:,i).^2)*tip_speed)/sound_speed);
end
inv_stall = zeros(length(r_segn),length(psi)); % inverse flow region and stalled stations
alpha_eff = NaN(length(r_segn),length(psi));   % non-stalled region of the rotor blade.
stall = NaN(length(r_segn),length(psi));       % positive stall region
for i = 1:length(r_segn)
    for j = 1:length(psi)
        if (alpha_e(i,j) - convang(alpha_stall_up,'deg','rad') > 0 && ...
                alpha_e(i,j) - convang(alpha_stall_up,'deg','rad') < 1)
            stall(i,j) = alpha_e(i,j);
        elseif (alpha_e(i,j) - convang(alpha_stall_up,'deg','rad') > 0 || ...
                alpha_e(i,j) - convang(alpha_stall_lo,'deg','rad') < 0)
            inv_stall(i,j) = alpha_e(i,j);
        else
            alpha_eff(i,j) = alpha_e(i,j);
        end
    end
end

% Compute first stalled sections and shows on screen warning
[ir,ic] = find((alpha_e(:,:) - convang(alpha_stall_up,'deg','rad') > 0) & ...
    (alpha_e(:,:) - convang(alpha_stall_up,'deg','rad') < convang(0.5,'deg','rad')));
if (~isempty(ir))
    warning('Rotor is stalling!')
end
% Compute max values for angles and Mach number
alpha_max = round(convang(max(alpha_eff,[],'all'),'rad','deg'),2);
Meff_max = round(max(Meff,[],'all'),2);

% Inizializing graphics' variables
mu = mu/cos(convang(alpha_tpp,'deg','rad')); % Dimensionless flight speed
% Creation of the mesh grid
[r_2d,psi_2d] = meshgrid(linspace(0,r_segn(1),length(r_segn)),psi);
x_hub = r_2d.*cos(psi_2d-pi/2);
y_hub = r_2d.*sin(psi_2d-pi/2);
[r_2d,psi_2d] = meshgrid(r_segn,psi);
x = r_2d.*cos(psi_2d);
y = r_2d.*sin(psi_2d);

%% Graphics

switch color
    case 'on'
        ifig = 1; % Indexing for figure
        % Mach colormap
        figure(ifig)
        set(figure(ifig),'Color','w');
        clf
        h = polar(x,y);
        hold on;
        fill(x_hub,y_hub,'k'); hold on; axis equal; axis off; grid off
        str = {'M_{e}(r,\psi)',['\mu = ' num2str(mu)], ...
            ['M_{e_{max}} = ' num2str(Meff_max)]};
        text(r_segn(end-10),r_segn(end-2),str,'Color','k','FontSize',10);
        contourf(x,y,Meff',40);
        axis equal; grid on; colormap jet;
        cb = colorbar('westoutside');
        position = cb.Position;
        cb.Position(1) = 0.55*position(1);
        cb.Position(2) = 1.7*position(2);
        cb.Position(3) = 0.5*position(3);
        cb.Position(4) = 0.8*position(4);
        cb.Label.String = ('M_e');
        cb.Label.FontSize = (10);
        set(h,'Visible','off')
        axis off; axis image
        view([90 90])
        ifig = ifig + 1;
        if (isempty(ir))
            % Alpha colormap
            figure(ifig)
            set(figure(ifig),'Color','w');
            clf
            h = polar(x,y);
            hold on;
            fill(x_hub,y_hub,'w'); hold on; axis equal; axis off; grid off
            str = {'\alpha_{e}(r,\psi)',['\mu = ' num2str(mu)], ...
                ['\alpha_{e_{max}} = ' num2str(alpha_max) '°']};
            text(r_segn(end-10),r_segn(end-2),str,'Color','k','FontSize',10);
            contourf(x,y,convang(alpha_eff','rad','deg'),40);
            axis equal; grid on;
            colormap jet;
            cb = colorbar('westoutside');
            position = cb.Position;
            cb.Position(1) = 0.55*position(1);
            cb.Position(2) = 1.7*position(2);
            cb.Position(3) = 0.5*position(3);
            cb.Position(4) = 0.8*position(4);
            cb.Label.String = ('\alpha_e(°)');
            cb.Label.FontSize = (10);
            set(h,'Visible','off')
            axis off; axis image
            view([90 90])
        else
            % Stalled sections with alpha color contour
            figure(ifig)
            set(figure(ifig),'Color','w');
            clf
            polar(psi(ic),r_segn(ir),'k*');
            hold on
            fill(x_hub,y_hub,'k'); hold on; axis equal; axis off; grid off
            str = {'\alpha_{e}(r,\psi)','* Stall',['\mu = ' num2str(mu)]};
            text(r_segn(end-10),r_segn(end-2),str,'Color','k','FontSize',10);
            contourf(x,y,convang(alpha_eff','rad','deg'),40);
            axis equal; grid on; axis off; axis image; colormap jet;
            cb = colorbar('westoutside');
            position = cb.Position;
            cb.Position(1) = 0.55*position(1);
            cb.Position(2) = 1.7*position(2);
            cb.Position(3) = 0.5*position(3);
            cb.Position(4) = 0.8*position(4);
            cb.Label.String = ('\alpha_e(°)');
            cb.Label.FontSize = (10);
            view([90 90])
        end
    case 'off'
        ifig = 1; % Indexing for figure
        % Mach Map
        figure(ifig)
        set(figure(ifig),'Color','w');
        clf
        h = polar(x,y);
        hold on;
        fill(x_hub,y_hub,'k'); hold on; axis equal; axis off; grid off
        str = {'M_{e}(r,\psi)',['\mu = ' num2str(mu)], ...
            ['M_{e_{max}} = ' num2str(Meff_max)]};
        text(r_segn(end-10),r_segn(end-2),str,'Color','k','FontSize',10);
        contour(x,y,Meff','k','ShowText','on');
        axis equal; grid on;
        set(h,'Visible','off')
        axis off; axis image
        view([90 90])
        ifig = ifig + 1;
       if (isempty(ir))
            % Alpha Map
            figure(ifig)
            set(figure(ifig),'Color','w');
            clf
            h = polar(x,y);
            hold on;
            fill(x_hub,y_hub,'k'); hold on; axis equal; axis off; grid off
            str = {'\alpha_{e}(r,\psi)',['\mu = ' num2str(mu)], ...
                ['\alpha_{e_{max}} = ' num2str(alpha_max) '°']};
            text(r_segn(end-10),r_segn(end-2),str,'Color','k','FontSize',10);
            contour(x,y,convang(alpha_eff','rad','deg'),'k','ShowText','on');
            axis equal; grid on;
            set(h,'Visible','off')
            axis off; axis image
            view([90 90])
            ifig = ifig + 1;
        else
            % Stall sections with alpha isolines
            figure(ifig)
            set(figure(ifig),'Color','w');
            clf
            polar(psi(ic),r_segn(ir),'k*');
            hold on
            fill(x_hub,y_hub,'k'); hold on; axis equal; axis off; grid off
            str = {'\alpha_{e}(r,\psi)','* Stall',['\mu = ' num2str(mu)]};
            text(r_segn(end-10),r_segn(end-2),str,'Color','k','FontSize',10);
            contour(x,y,convang(alpha_eff','rad','deg'),'k','ShowText','on');
            axis equal; grid on; axis off; axis image;
            view([90 90])
        end
end
end
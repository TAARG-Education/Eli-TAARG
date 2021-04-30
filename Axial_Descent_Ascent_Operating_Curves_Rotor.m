%% \Axial_Descent_Ascent_Operating_Curves_Rotor.m
%  \brief: the function plots w(V_infty) and P(V_infty) curves according to
%  Impulsive theory.
%  aerodynamic model
%  \author: Colledà Moreno, Veneruso Salvatore
%  \version: 1.00
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
% |Name        : Axial_Descent_Ascent_Operating_Curves_Rotor.m                           |
% |Author      : Colledà Moreno, Veneruso Salvatore                                      |
% |              University of Naples Federico II.                                       |
% |Version     : 1.01                                                                    |
% |Date        : 21/12/2020                                                              |
% |Modified    : 08/01/2021                                                              |
% |Description : the function plots w(V_infty) and P(V_infty) curves according to
%                Rotor Simply Impulsive theory.                                          |
% |Reference   : Renato Tognaccini. Appunti Aerodinamica dell'ala rotante.
%                Università degli studi di Napoli Federico II. a.a.2020/2021             |
% |Input       : * the inputs must be Mass of rotorcraft and radius of rotor             |
% |Output      : w(V_infty) and P(V_infty) plots                                         |
% |Note        :                                                                         |
% ==============================================================================================
function [Power,Induction] =Axial_Descent_Ascent_Operating_Curves_Rotor(M,R,hnew,V_inf)

g = 9.81;
%Altitude Range
hh    = linspace(0,6000,100);

%Axial Induction

[~, ~, ~, rho1] = atmosisa(hh);

rho = @(h) interp1(hh,rho1,h, 'pchip');
wh  = @(h) sqrt((M*g)/(2*rho(h)*pi*R^2));

% Non-Dimensional Variables Definition
V_tilde = @(V,h) V/wh(h);

w_tilde_salita  = @(V,h) (-V_tilde(V,h)/2) + sqrt((V_tilde(V,h)/2)^2+1);
w_tilde_discesa = @(V,h) (-V_tilde(V,h)/2) - sqrt((V_tilde(V,h)/2)^2-1);

if nargin==3
    VVs   = linspace(0, 50, 100);
    VVd   = linspace(-50, 0, 7000);
    
    %Non Dimensional Induction
    
    WTS = zeros(length(hh),length(VVs)); %Ascent Induction
    WTD = zeros(length(hh),length(VVd)); %Descent Induction
    aa  = zeros(1,length(hh));           %Control Parameter
    
    % Matrices fill
    for i = 1 : length(hh)
        for j = 1 : length(VVs)
            WTS(i,j) = w_tilde_salita(VVs(j),hh(i));
        end
    end
    
    for i = 1 : length(hh)
        for j = 1 : length(VVd)
            if V_tilde(VVd(j),hh(i)) <= -2            %Validity limit of simply impulsive theory;
                WTD(i,j) = w_tilde_discesa(VVd(j),hh(i));
                aa(1,i) = j;
            else
                WTD(i,j) = 0;
            end
        end
    end
    
    %Non-Dimensional Power
    P_tilde_salita   = @(V,h) V_tilde(V,h) + w_tilde_salita(V,h);
    P_tilde_discesa  = @(V,h) V_tilde(V,h) + w_tilde_discesa(V,h);
    PTS = zeros(length(hh),length(VVs)); %Ascent Power
    PTD = zeros(length(hh),length(VVd)); %Descent Power
    bb  = zeros(1,length(hh));           %Control Paramenter
    
    % Matrices fill
    for i = 1 : length(hh)
        for j = 1 : length(VVs)
            PTS(i,j) = P_tilde_salita(VVs(j),hh(i));
        end
    end
    
    for i = 1 : length(hh)
        for j = 1 : length(VVd)
            if V_tilde(VVd(j),hh(i)) <= -2            %Validity limit of simply impulsive theory;
                PTD(i,j) = P_tilde_discesa(VVd(j),hh(i));
                bb(1,i) = j;
            else
                break
            end
        end
    end
    
    %% 3D Plots of Induction [output]
    figure(1)
    plot3(hh(1)*ones(1,length(WTS(1,:))),VVs,WTS(1,:),'k')
    hold on
    for i = 2 : length(hh)
        plot3(hh(i)*ones(1,length(WTS(i,:))),VVs,WTS(i,:),'k')
    end
    for i = 1 : length(hh)
        plot3(hh(i)*ones(1,length(WTD(i,1:aa(i)))),VVd(1,1:aa(i)),WTD(i,1:aa(i)),'k')
    end
    yticks([-50 -25 0 25 50])
    xlabel('$Height[m]$','Interpreter','latex', 'FontSize',15)
    ylabel('$V_{\infty}[\frac{m}{s}]$','Interpreter','latex', 'FontSize',15)
    zlabel('$\widetilde{w}$','Interpreter','latex', 'FontSize',15)
    grid on
    title('$w (V_{\infty}, h)$','Interpreter','latex', 'FontSize', 15)
    view(71,32)
    
    % 3D Plots of Power [output]
    figure(2)
    plot3(hh(1)*ones(1,length(PTS(1,:))),VVs,PTS(1,:),'k')
    hold on
    for i = 2 : length(hh)
        plot3(hh(i)*ones(1,length(PTS(i,:))),VVs,PTS(i,:),'k')
    end
    for i = 1 : length(hh)
        plot3(hh(i)*ones(1,length(PTD(i,1:bb(i)))),VVd(1,1:bb(i)),PTD(i,1:bb(i)),'k')
    end
    yticks([-50 -25 0 25 50])
    xlabel('$Height[m]$','Interpreter','latex', 'FontSize',15)
    ylabel('$V_{\infty}[\frac{m}{s}]$','Interpreter','latex', 'FontSize',15)
    zlabel('$\widetilde{P}$','Interpreter','latex', 'FontSize',15)
    grid on
    title('$P (V_{\infty}, h)$','Interpreter','latex', 'FontSize', 15)
    view(71,32)
    
    %% Insertion of Interes Altitude;
    
    WTSnew = zeros(1,length(VVs));  %Vector Inizialization of axial ascent induction related to velocity
    WTDnew = zeros(1,length(VVd));  %Vector Inizialization of axial descent induction related to velocity
    
    % Matrices fill
    for j = 1 : length(VVs)
        WTSnew(1,j) = w_tilde_salita(VVs(j),hnew);
    end
    
    for j = 1 : length(VVd)
        if V_tilde(VVd(j),hnew) <= -2                  %Validity limit of simply impulsive theory;
            WTDnew(1,j) = w_tilde_discesa(VVd(j),hnew);
            aanew = j;
        else
            break
        end
    end
    
    PTSnew = zeros(1,length(VVs));
    PTDnew = zeros(1,length(VVd));
    
    
    for j = 1 : length(VVs)
        PTSnew(1,j) =  P_tilde_salita(VVs(j),hnew);
    end
    
    for j = 1 : length(VVd)
        if V_tilde(VVd(j),hnew) <= -2                  %Validity limit of simply impulsive theory;
            PTDnew(1,j) = P_tilde_discesa(VVd(j),hnew);
            bbnew = j;
        else
            break
        end
    end
    
    %% 2D Plots of induction [Output]
    figure(3)
    plot(V_tilde(VVs,hnew), WTSnew, '-k')
    hold on
    plot(V_tilde(VVd(1:aanew),hnew), WTDnew(1:aanew), '-k')
    grid on
    xlabel('$\widetilde{V}$','Interpreter','latex','FontSize',15)
    ylabel('$\widetilde{w}$','Interpreter','latex','FontSize',15)
    title(['$\widetilde{w}$ per h=' num2str(hnew) '[m]'],'Interpreter','latex')
    
    % 2D Plots of power [Output]
    figure(4)
    plot(V_tilde(VVs,hnew), PTSnew, '-k')
    hold on
    plot(V_tilde(VVd(1:aanew),hnew), PTDnew(1:bbnew), '-k')
    grid on
    xlabel('$\widetilde{V}$','Interpreter','latex','FontSize',15)
    ylabel('$\widetilde{P}$','Interpreter','latex','FontSize',15)
    title(['$\widetilde{P}$ per h=' num2str(hnew) '[m]'],'Interpreter','latex')
    %% Numerical Output of Power and Induction for the interest altitude;
    Power = [PTDnew(1:bbnew), PTSnew];
    Induction = [WTDnew(1:aanew), WTSnew];
    
elseif nargin==4 % In Output only value of Power and Induction
    % at altitude and velocity of interest;
    
    if V_tilde(V_inf,hnew) >= 0
        w = w_tilde_salita(V_inf,hnew)*wh(hnew);
        P   = (M*g)*(V_inf + w);
    elseif V_tilde(V_inf,hnew)<= -2
        w= w_tilde_discesa(V_inf,hnew)*wh(hnew);
        P  = (M*g)*(V_inf + w);
    else
        %output
        errordlg('Simply Impulsive Rotor Theory not respected','ERROR!!');
        %Output
        Power = [];
        Induction = [];
        return
    end
    %output
    ff = msgbox(sprintf('Power= %d [kW], \n Induction= %d [m/s]', P/1000, w),...
        'Power and Induction at the altitude and velocity of interest');
    set(ff, 'position', [500 250 400 65]);
    %Output
    Power = P;
    Induction = w;
end

end
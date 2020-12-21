function [Power,Induction] = Axial_Descent_Ascent_Operating_Curves_Rotor(M,R)
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
% |Version     : 1.00                                                                    |
% |Date        : 21/12/2020                                                              |
% |Modified    : 21/12/2020                                                              |
% |Description : the function plots w(V_infty) and P(V_infty) curves according to
%                Rotor Simply Impulsive theory.                                          |     
% |Reference   : Renato Tognaccini. Appunti Aerodinamica dell'ala rotante.
%                Università degli studi di Napoli Federico II. a.a.2020/2021             |
% |Input       : * the inputs must be Mass of rotorcraft and radius of rotor             | 
% |Output      : w(V_infty) and P(V_infty) plots                                         |
% |Note        :                                                                         |
% ==============================================================================================

%% Function Main

%Definizione Variabili
hh    = linspace(0,6000,1000);
VVs   = linspace(0, 50, 100);
VVd   = linspace(-50, 0, 100);
rho1  = zeros(1,length(hh));
g     = 9.81;
% wh  = zeros(1,length(hh));

%Calcolo Induzione in hovering

[~, ~, ~, rho1] = atmosisa(hh);
    
rho = @(h) interp1(hh,rho1,h, 'pchip');
wh  = @(h) sqrt((M*g)/(2*rho(h)*pi*R^2));

% Definizione variabili adimensionali 
V_tilde = @(V,h) V/wh(h);

%Calcolo curve di Funzionamento:

%Induzione
w_tilde_salita  = @(V,h) (-V_tilde(V,h)/2) + sqrt((V_tilde(V,h)/2)^2+1);
w_tilde_discesa = @(V,h) (-V_tilde(V,h)/2) - sqrt((V_tilde(V,h)/2)^2-1);

WTS = zeros(length(hh),length(VVs)); %Induzione in salita
WTD = zeros(length(hh),length(VVd)); %Induzione in discesa
aa  = zeros(1,length(hh));


for i = 1 : length(hh)
    for j = 1 : length(VVs)
        WTS(i,j) = w_tilde_salita(VVs(j),hh(i));
    end
end

for i = 1 : length(hh)
    for j = 1 : length(VVd)
        if V_tilde(VVd(j),hh(i)) < -2
        WTD(i,j) = w_tilde_discesa(VVd(j),hh(i));
        aa(1,i) = j;
        else 
            WTD(i,j) = 0;
        end
    end
end

%Potenza

P_tilde_salita   = @(V,h) V_tilde(V,h) + w_tilde_salita(V,h);
P_tilde_discesa  = @(V,h) V_tilde(V,h) + w_tilde_discesa(V,h);

PTS = zeros(length(hh),length(VVs)); %Induzione in salita
PTD = zeros(length(hh),length(VVd)); %Induzione in discesa
bb  = zeros(1,length(hh));

for i = 1 : length(hh)
    for j = 1 : length(VVs)
        PTS(i,j) = P_tilde_salita(VVs(j),hh(i));
    end
end

for i = 1 : length(hh)
    for j = 1 : length(VVd)
        if V_tilde(VVd(j),hh(i)) < -2
        PTD(i,j) = P_tilde_discesa(VVd(j),hh(i));
        bb(1,i) = j;
        else 
            break
        end
    end
end


%%
% Plot Curva funzionamento Induzione funzione di velocità e quota
figure(1)
 plot3(hh(1)*ones(1,length(WTS(1,:))),VVs,WTS(1,:))
 hold on
 for i = 2 : length(hh)
     plot3(hh(i)*ones(1,length(WTS(i,:))),VVs,WTS(i,:))
 end
 for i = 1 : length(hh)
     plot3(hh(i)*ones(1,length(WTD(i,1:aa(i)))),VVd(1,1:aa(i)),WTD(i,1:aa(i)))
 end
xlabel('$Quota$','Interpreter','latex', 'FontSize',15)
ylabel('$V_{\infty}$','Interpreter','latex', 'FontSize',15)
zlabel('$w$','Interpreter','latex', 'FontSize',15)
grid on
title('$w (V_{\infty}, h)$','Interpreter','latex', 'FontSize', 15)
view(71,32)

% Plot Curva funzionamento Potenza funzione di velocità e quota
figure(2)
 plot3(hh(1)*ones(1,length(PTS(1,:))),VVs,PTS(1,:))
 hold on
 for i = 2 : length(hh)
     plot3(hh(i)*ones(1,length(PTS(i,:))),VVs,PTS(i,:))
 end
 for i = 1 : length(hh)
     plot3(hh(i)*ones(1,length(PTD(i,1:bb(i)))),VVd(1,1:bb(i)),PTD(i,1:bb(i)))
 end
xlabel('$Quota$','Interpreter','latex', 'FontSize',15)
ylabel('$V_{\infty}$','Interpreter','latex', 'FontSize',15)
zlabel('$P$','Interpreter','latex', 'FontSize',15)
grid on
title('$P (V_{\infty}, h)$','Interpreter','latex', 'FontSize', 15)
view(71,32)

%% Scelta Quota di Interesse 
prompt = {'Inserire la quota di interesse in metri [min=0,Max=6000]: '};
dlgtitle = 'Quota';
dims = [1 35];
answer = inputdlg(prompt,dlgtitle,dims);
hnew = str2double(answer{1});
% disp('Premere un tasto per scegliere la quota di interesse...')
% pause
% hnew = input('Inserire la quota di interesse: ');
WTSnew = zeros(1,length(VVs));
WTDnew = zeros(1,length(VVd));


for j = 1 : length(VVs)
    WTSnew(1,j) = w_tilde_salita(VVs(j),hnew);
end

for j = 1 : length(VVd)
    if V_tilde(VVd(j),hh(i)) < -2
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
    if V_tilde(VVd(j),hh(i)) < -2
        PTDnew(1,j) = P_tilde_discesa(VVd(j),hnew);
        bbnew = j;
    else
        break
    end
end

%% Plot Curve A Quota Di Interesse
figure(3)
plot(VVs, WTSnew, '-k')
hold on
plot(VVd(1:aanew), WTDnew(1:aanew), '-k')
grid on
xlabel('$V_{\infty}$','Interpreter','latex','FontSize',15)
ylabel('$w$','Interpreter','latex','FontSize',15)
title(['$w$ per h=' num2str(hnew)],'Interpreter','latex')


figure(4)
plot(VVs, PTSnew, '-k')
hold on
plot(VVd(1:aanew), PTDnew(1:bbnew), '-k')
grid on
xlabel('$V_{\infty}$','Interpreter','latex','FontSize',15)
ylabel('$P$','Interpreter','latex','FontSize',15)
title(['$P$ per h=' num2str(hnew)],'Interpreter','latex')

%Per ora li ho messi a caso per far funzionare il tutto.

Power = 1;
Induction = WTS;




end
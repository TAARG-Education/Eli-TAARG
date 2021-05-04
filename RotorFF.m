%% \RotorFF
%  \The function plots the characteristic curves for rotor in forward
%  flight for both constant thrust and power and gives in output also the
%  relative x and y values. It requires in input the angle of attack in 
%  degrees.
%  The plot avaiable are:
%- For constant Thrust:
%  -w versus V
%  -P versus V
%- For constant Power:
%  -w versus V
%  -T versus V
% where w = induction, V = asymptotic velocity, T = Thrust, P = Power
% All the values are non-dimensional with respect to their value in
% hovering (for V is used induction in hovering).
%  \Angelo Tagliafierro
%  \0.1
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
% |Name        : RotorFF                                                                       |
% |Author      : Angelo Tagliafierro                                                           |
% |              University of Naples Federico II.                                             |
% |Version     : 0.1                                                                           |
% |Date        : 15/12/2020                                                                     |
% |Modified    : 15/12/2020                                                                    |
% |Description : The function generates the characteristic curves for rotor in forward         |
%                flight and the relative x and y values                                        |
% |Reference   : R. Tognaccini. "Lezioni di aerodinamica dell’ala rotante" 2019 pp. 84-85.     |
% |Input       : Angle of attack                                                               |
% |Output      : Plots and axis values                                                         |
% |Note        :                                                                               |
% ==============================================================================================

function [Vt,wt,Vp,wp,Pt,Tp] = RotorFF(alfa)

%%Angle conversion to radiant
alfa=deg2rad(alfa);

%% Constant thrust
%Function definition
f1 = @(x,y) (((x.*cos(alfa)).^2 + (x.*sin(alfa)+y).^2).*(y.^2))-1;

%Plot setting
subplot(2,2,1);
h1=fimplicit(f1,[0 10 0 1],'k');
axis square;
xlabel('$\tilde{V_{\infty}}$','Interpreter','latex');ylabel('$\tilde{w}$','Interpreter','latex');
title('Characteristic curve $\tilde{w} = \tilde{w}(\tilde{V_{\infty}})$,constant Thrust', 'Interpreter','latex') ;

%Axis values 
Vt = (get(h1, 'XData'));
wt = (get(h1, 'YData'));

%Function of Power versus speed for costant Thrust
Pt=Vt.*sin(alfa)+wt;

%Plot setting
subplot(2,2,3);
plot(Vt,Pt,'k')
axis([0 max(Vt) 0 max(Pt)]);
axis square;
xlabel('$\tilde{V_{\infty}}$','Interpreter','latex');ylabel('$\tilde{P}$','Interpreter','latex');
title('Characteristic curve $\tilde{P} = \tilde{P}(\tilde{V_{\infty}})$,constant Thrust', 'Interpreter','latex') ;


%% Constant Power

%Function Definition
f2 = @(x,y) (((x.*cos(alfa)).^2 + (x.*sin(alfa)+y).^2)).*(y.^2).*((x.*sin(alfa)+y))-1 ;

%Plot setting
subplot(2,2,2);
h2=fimplicit(f2,[0 5 0 1],'k');
axis square;
xlabel('$\tilde{V_{\infty}}$','Interpreter','latex');ylabel('$\tilde{w}$','Interpreter','latex');
title('Characteristic curve $\tilde{w} = \tilde{w}(\tilde{V_{\infty}})$,constant Power', 'Interpreter','latex') ;

%Axis values
Vp= (get(h2, 'XData'));
wp= (get(h2, 'YData'));

%Function of Thrust versus speed for costant Power
Tp=(Vp*sin(alfa)+wp).^(-1);

%Plot setting
subplot(2,2,4);
plot(Vp,Tp,'k')
axis([0 max(Vp) 0 max(Tp)]);
axis square;
xlabel('$\tilde{V_{\infty}}$','Interpreter','latex');ylabel('$\tilde{T}$','Interpreter','latex');
title('Characteristic curve $\tilde{T} = \tilde{T}(\tilde{V_{\infty}})$,constant Power', 'Interpreter','latex') ;

end

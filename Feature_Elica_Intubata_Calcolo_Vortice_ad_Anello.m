function [fx] = RVortexInt(ch,Di,De)
%% \RVortexInt.m
%  \The function evaluates the dimensionless radial velocity component
%   of the ring vortex by means of integration of Biot-Savar law
%  \Olino Massimiliano, Marino Giuseppe
%  \1.0
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
% |Name        : RVortexInt.m                                                                  |
% |Author      : Olino Massimiliano, Marino Giuseppe                                           |
% |              University of Naples Federico II.                                             |
% |Version     : 1.0                                                                           |
% |Date        : data di creazione                                                             |
% |Modified    : data ultima modifica                                                          |
% |Description : The function evaluates the dimensionless radial velocity component            |
% |              of the ring vortex by means of integration of Biot-Savar law                  |
% |Reference   : McCormick, B.W.,(1967), Aerodynamics of V/STOL Flight, Academic Press.        |
% |              Tognaccini, R., (2020), Lezioni di AERODINAMICA DELL'ALA ROTANTE.             |
% |Input       : ch (chord), Di (inlet diameter), De (exit diameter) of the shroud             |
% |Output      : fx (the dimensionless radial velocity component of the ring vortex )          |
% |Note        :                                                                               |
% ==============================================================================================
%
%Gui for the input by users.
title='Input Data';
prompt={'Shroud chord [m]','Inlet diameter [m]','Exit diameter [m]'};
answer=inputdlg(prompt,title,[1 55;1 55;1 55]);
ch=str2double(answer{1});
Di=str2double(answer{2});
De=str2double(answer{3});

% Creation of chroud geometry
t=linspace(0,1,1000);      
    for i=1:length(t) 
        xs(i)=0+t(i)*(ch);
        ys(i)=Di/2+t(i)*(De/2-Di/2);
    end
D14=2*ys(250);               % Diameter at 1/4 of the shroud
D34=2*ys(750);               % Diameter at 3/4 of the shroud
c_D14=ch/D14 ;
RatioD=D34/D14;

r=(D14/2);                   % flow radius
x=D34/2;
z=ch*(3/4);
y=0;

R1=@(t) cos(t)./(((x-r.*cos(t)).^2 +((y-r.*sin(t)).^2)+z^2)).^(3/2);
R2=@(t) sin(t)./(((x-r.*cos(t)).^2 +(y-r.*sin(t)).^2+z^2)).^(3/2);
R3=@(t) ((x-r.*cos(t)).*cos(t)+(y-r.*sin(t)).*sin(t))./(((x-r.*cos(t)).^2 +(y-r.*sin(t)).^2+z^2)).^(3/2);
Ix=integral(R1,0,2*pi);
Iy=integral(R2,0,2*pi);
Iz=integral(R3,0,2*pi);

% Adimesionalisation of the velocity components

fx=(r*z*D14/4)*Ix;           % Component of interests
fy=-(r*z*D14/4)*Iy;
fz=-(r*D14/4)*Iz;

% Vector sum
f=sqrt(fx^2+fy^2+fz^2);

% Results
figure(1);
plot(f,f)
axis ([0 1 0 1]);
text(0.25,1,'Velocity induced by vortex ring');
text(0.25,0.90,['fx= ',num2str(fx)]);
axis off;
end


%% \file autogiro_performances
%  \brief Autogiro performances
%  \author Filippo Borriello - Alessandro Furno
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
% |Name        : autogiro_performances                                                         |
% |Author      : Filippo Borriello - Alessandro Furno                                          |
% |              University of Naples Federico II.                                             |
% |Version     : 1.0.0                                                                         |
% |Date        : 16/02/2023                                                                    |
% |Modified    : ----------                                                                    |
% |Description : This function provides forward flight performances                            |
%                of an autogiro: the procedure is based upon the one shown                     |
%                in reference. For a given autogiro, intake ratio is to                        |
%                be set by the user.                                                           |
%                Autorotation condition yields to equations (see                               |
%                reference, page 110) which are here dealt with                                |
%                by means of the Matlab Symbolic Toolbox. Once solved                          |
%                in respect to the the collective pitch, expressions of thrust                 |
%                coefficient, angle of attack, angular velocity,                               |
%                freestream velocity can be evaluated. Please note that                        |
%                it is not sure to find a solution: it simply means that                       |
%                autorotation cannot be achieved under the assigned                            |
%                conditions. This is not taken into account directly in                        |
%                this version, so the user should pay attention to give                        |
%                resonable value of the advance ratio.                                         |
%                                                                                              |
% |Reference   : Lezioni di Aerodinamica dell'Ala Rotante a.a. 2022/2023                       |
%                Renato Tognaccini                                                             |
% |Input       : Autogiro mass, W [kg]                                                         |
%                Rotorblade radius, R [m]                                                      |
%                Number of blades, N                                                           |
%                Chord length, c [m]                                                           |
%                Lock's number, gamma                                                          |
%                Blade twisting, theta_tw [rad] (assuming linear variation)                    |
%                Equivalent wet area, f [m^2]                                                  |
%                Advance ratio, mu                                                             |
%                Descent angle, X [deg]                                                        |
% |Output      : Thrust coefficient, T_C                                                       |
%                Rotor drag coefficient, H_C                                                   |
%                Angle of attack, alpha_deg [deg]                                              |
%                Angular velocity, omega [rad/s]                                               |
%                Freestream velocity,v0 [m/s]                                                  |
%                Collective pitch, Theta_0 [°]                                                 |
% |Note        : note addizionali                                                              |
% ==============================================================================================


function [pnec, alpha, omega,THETA_0F,v0, pi, rs ]= autogiro_performances(w,ro,co,cla,r,tetatw,n,fd,A,cd0,lock,lambda,mu,Pmax)
sigma=n*co/(3.14*r);
syms teta0 omega [1 50]
for i=1:50
    % thrust coefficient
Tc(i) = 0.5*sigma*cla*(teta0(i)/3*(1+3/2*mu(i)^2)+tetatw/4*(1+mu(i)^2)-lambda(i)/2);   

 % parasitic torque coefficient
Qc0(i) = sigma*cd0/8*(1+mu(i)^2);  

% flapping coefficients expressions
b_0(i) = lock*(teta0(i)/8*(1+mu(i)^2)+tetatw/10*(1+5/6*mu(i)^2)-lambda(i)/6);         
b_1c(i) = -2*mu(i)*(4/3*teta0(i)+tetatw-lambda(i))/(1-0.5*mu(i)^2);
b_1s(i) = -4/3*mu(i)*b_0(i)/(1+0.5*mu(i)^2);

 % induced rotor drag coefficient
H_ci(i) = 0.5*sigma*cla*(teta0(i)*(-1/3*b_1c(i)+0.5*mu(i)*lambda(i))...                
+tetatw*(-1/4*b_1c(i)+1/4*mu(i)*lambda(i))+3/4*lambda(i)*b_1c(i)...
+1/6*b_0(i)*b_1s(i)+1/4*mu(i)*(b_0(i)^2+b_1c(i)^2));

% parasitic rotor drag coefficient
H_C0(i) = sigma*cd0*mu(i)/4;     

%from the definitiion of autorotaion
Pc(i) = lambda(i)*Tc(i)+Qc0(i)-mu(i)*H_ci(i);                                           

% we choose the second solotuion
THETA_0(i,:) = vpasolve(Pc(i)==0,teta0(i));

THETA_0F(i)=THETA_0(i,2);

% It is possible now to find Tc=Tc(mu),Hc=Hc(mu),alfa=alfa(mu)

TC(i) = 0.5*sigma*cla*(THETA_0F(i)/3*(1+3/2*mu(i)^2)...
+tetatw/4*(1+mu(i)^2) -lambda(i)/2);

B_0(i) = lock*(THETA_0F(i)/8*(1+mu(i)^2)+tetatw/10*(1+5/6*mu(i)^2)...
-lambda(i)/6);

B_1c(i)= -2*mu(i)*(4/3*THETA_0F(i)+tetatw-lambda(i))/(1-0.5*mu(i)^2);

B_1s_lambda(i) = -4/3*mu(i)*B_0(i)/(1+0.5*mu(i)^2);

H_Ci(i)= 0.5*sigma*cla*(THETA_0F(i)*(-1/3*B_1c(i)...
+0.5*mu(i)*lambda(i))+tetatw*(-1/4*B_1c(i)+1/4*mu(i)*lambda(i))...
+3/4*lambda(i)*B_1c(i)+1/6*B_0(i)*B_1s_lambda(i)...
+1/4*mu(i)*(B_0(i)^2+B_1c(i)^2));

 % induced intake ratio
lambda_i(i) = 1.2*TC(i)/(2*sqrt(mu(i)^2+lambda(i)^2));                                    

 % angle of attack
alpha(i) = atan(1/mu(i)*(lambda(i)-TC(i)/(2*sqrt(mu(i)^2+lambda(i)^2))));                 

 % angular velocity of the rotor
omega(i) =sqrt(w/(cos(alpha(i))*TC(i)*ro*3.14*r^4));                                      

end

%Evaluation of the output

for i=1:50
lambda_i(i)=vpa(lambda_i(i));

 % Thrust coefficient
TC(i)=vpa(TC(i));                                            

% [rad/s] Angular velocity
omega(i)=vpa(omega(i));                                       

% Rotor drag coefficient                                              
H_Ci(i)=vpa(H_Ci(i));

H_C0(i)=vpa(H_C0(i));

HC(i)=(H_Ci(i)+H_C0(i));

% Collective pitch
THETA_0F(i)=vpa(THETA_0F(i))/0.0174;                          

alpha(i)=vpa(alpha(i))/0.0174; 

% freestream velocity
v0(i)=mu(i)*omega(i)*r;                                       

 % induced power coefficient
pi(i)=TC(i)*lambda_i(i)*ro*omega(i)^3*r^5*3.14;              

% fuselage parastic coefficient
pfus(i)=0.5*fd*v0(i)^3*A;                                     

prot(i)=mu(i)*(H_Ci(i)+H_C0(i))*ro*3.14*r^5*omega(i)^3;

% total required power
pnec(i)=((pi(i)+prot(i))+pfus(i))/0.7;                        

vprot(i)=mu(i)*tan(alpha(i)*0.01744);

% climb rate
rs(i)=(Pmax-pnec(i))/w*3.28*60;                               

end
end

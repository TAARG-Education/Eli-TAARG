clear; clc; close all;

% This is the main code written to validate the function.

%Input data
        Cla=6.28;           %
        alphamax = 14;      %Maximun angle of attack
        c = 1;              %Chord
        R = 30;             %Radius
        N = 3;              %Number of blade
        cd = 0.01;          %Drag coefficient
            
 
 load DeVries_data.mat %import data from De Vries
 
 %Function 
 [Cp, Lambda, Cq] = Turbine_single_tube(Cla, c, R, N, cd, alphamax);
 
        
%Ideal condition
lambdai=linspace(0.01,12,100);
ai=N.*c.*Cla.*lambdai./(16.*R);
cpi=4.*ai.*(1-ai).^2-N.*c.*cd.*lambdai.^3./(2.*R);

%Plot cp, lambda
figure (1);
plot(lambda1,cp1,'kv'); hold on; grid on; %Devries data
plot(Lambda,Cp,'k.',lambdai,cpi,'k--');
legend('Single streamtube theory deVrei','Single streamtube theory funzione','Asintotic Analisys','Location','best');
        h=xlabel('$\lambda$','Interpreter','latex'); set(h,'Fontsize',20);
        h=ylabel('$C_P$','Interpreter','latex');
        set(h,'Fontsize',20);
        title('\sigma = 0.1 - Cd = 0.01');
        
figure (2);
plot(Lambda,Cq,'k .');
legend('Single streamtube theory Cq','Location','best');
        axis on; grid on;
        h=xlabel('$\lambda$','Interpreter','latex'); set(h,'Fontsize',20);
        h=ylabel('$C_Q$','Interpreter','latex');
        set(h,'Fontsize',20);
        title('\sigma = 0.1 - Cd = 0.01');

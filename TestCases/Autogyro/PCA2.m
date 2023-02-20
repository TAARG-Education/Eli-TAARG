clear all; close all; clc;
% Test case
%% Pitcairn PCA-2 autogiro in forward flight

% input 

W = 1563*9.81;                % MTOW [ N ]
Pd = 250000;                  % Maximum aviable power [ W ]
N = 4;                        % blade numbers
R = 8;                     % rotorblade radius [m]
A = 3.14*R^2;                 % rotor area [m^2]
c=0.63;                       % chord length [m]
Cl_a=5.4;
tetatw=-10*0.01744;           % blade twisting (linear variation is assumed) [°]
fd=0.01;                      % equivalent wet area [m^2] 
rho = 1.225;                  % air density [kg/m^3]
cd0=0.010;                    % avarage drag coefficient
lock=8;                       % Lock's number

lambda(1:50)=-0.01;
mu=linspace(0.03,0.5,50);

%% Function output
[Pnec, alpha, omega,THETA_0F,v0,pi,rs ] = AUTOGIRO_FUN(W,rho,c,Cl_a,R,tetatw,N,fd,A,cd0,lock,lambda,mu,Pd);

P=Pd*(ones(1,50)); % for grapchis reason only

% figure(1)
% plot(mu,Pnec,'k','linewidth',1.5); grid on;
% xlabel('\mu'); ylabel('P_n');
% title('Required power')
% 
% figure(2)
% plot(v0,THETA_0F,'k','linewidth',1.5); grid on;
% xlabel('V0 m/s'); ylabel(' θ_0 °');
% title('Collective pitch \theta °')
% 
% figure(3)
% plot(v0,alpha,'k','linewidth',1.5); grid on;
% xlabel('V0 m/s'); ylabel('alpha °')
% axis([0 60 -90 0])
% title('Angle of attack \alpha °');
% 
figure(4)
plot(v0,Pnec,'k','linewidth',1.5); hold on;
plot(v0,P,'--k','linewidth',1.5); 
grid on;
axis([0 60 0 30e4])
ylabel('Required Power P_n'); xlabel('V0 m/s')

% figure(5)
% plot(mu,pi,'k','linewidth',1.5); grid on;
% xlabel('\mu');ylabel('P_i');title('Potenza Indotta');
% 
figure(6)
plot(v0,rs,'k','linewidth',1.5);grid on;
ylabel('Climb rate  ft/min')
xlabel('V0 m/s')
axis ([0 60 0 2500])
title('Climb rate')

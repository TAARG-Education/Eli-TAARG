function [T,TS] = elica_intubata(rho,Gamma,Vinf,Dquarter,R,TR)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


%Calcolo la velocità indotta sull'elica
%eq 4.24
w0 = .5*(-Vinf + sqrt(Vinf^2 + 2*TR/(rho*pi*R^2)));

%utilizzo l'eq della continuità per calcolare il raggio del tubo di flusso
%eq 4.20
A = pi*R^2;
const=0;
const = (Vinf + w0)*A;
rt = sqrt(const/(Vinf*pi));

%calcolo la velocità radiale indotta dal rotore a c/4
%eq 4.23
cquarter=1/4;
viRquarter = -.5*rt*w0*R^2/((R^2+cquarter^2)^1.5);

%calcolo il contributo di spinta fornito dalla carenatura
TS = -rho*viRquarter*Gamma*pi*Dquarter;

T = TR + TS;




end


clc;clear all; close all;

N=3;
test_geometry=load('testxrotor.txt');
r=test_geometry(:,1)*8.5;
c=test_geometry(:,2)*8.5;
beta=convang((test_geometry(:,3)),'deg','rad');


TEST_ROT=load('test.txt');


LAMBDA=1./TEST_ROT(:,2);
V=8.90;

[Cp, Cq, Ct, Lambda]=Curve_caratteristiche_aeromotori_asse_orizzontal(N,r,beta,c,LAMBDA(1),LAMBDA(end));


P=TEST_ROT(:,10)*1000;
T=TEST_ROT(:,11);
Q=TEST_ROT(:,12);

CP=-P/(0.5*1.23*V^3*pi*8.5^2);
CQ=-Q/(0.5*1.23*V^2*pi*8.5^3);
CT=-T/(0.5*1.23*V^2*pi*8.5^2);


figure(1)
plot(Lambda,Ct,'-k',LAMBDA,CT,'-or')
xlabel('\lambda'); ylabel('C_T');
grid on;
legend('Funzione','XRotor')
axis([LAMBDA(end) LAMBDA(1) 0.3 1])


figure(2)
plot(Lambda,Cq,'-k',LAMBDA,CQ,'-or')
xlabel('\lambda'); ylabel('C_Q');
grid on;
axis([LAMBDA(end) LAMBDA(1) 0 0.4])
legend('Funzione','XRotor')


figure(3)
plot(Lambda,Cp,'-k',LAMBDA,CP,'-or')
axis([LAMBDA(end) LAMBDA(1) 0 0.7])
xlabel('\lambda'); ylabel('C_P');
grid on;
legend('Funzione','XRotor')




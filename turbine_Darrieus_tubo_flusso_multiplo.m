clear vars; clc; close all;

%% Input
c = 1; %corda del profilo
R = 30; %raggio pala 
phiv = linspace(0,360,100); %dominio di phi (posizione angolare dell'ala)
phiv = deg2rad(phiv); % dominio di phi in rad
cla = 2*pi; %clalfa profilo
lambdav = linspace(2,15,100); %vettore dei lambda
N = 3; %numero delle pale
cd = 0.01; %coefficiente di resistenza medio del profilo
cp = zeros(numel(lambdav),1); %inizializzazione del vettore del coefficiente di potenza
cq = zeros(numel(lambdav),1); %inizializzazione del vettore del coefficiente di coppia
a = zeros(1,numel(phiv)); %inizializzazione del vettore di induzione
v_vinf = zeros(numel(phiv),numel(phiv)); %inizializzazione del vettore di v/vinfinito
alpha = zeros(numel(phiv),numel(phiv)); %inizializzazione del vettore degli alfa
ind = zeros(numel(phiv),numel(phiv));
%% Calcolo di a, CP, CQ
for j = 1 : numel(lambdav)
    lambda = lambdav(j); %assegnazione di lambda
for i = 1 : numel(phiv) 
    phi = phiv(i); %assegnazione di phi
    eq = @(a) ((1-a).*a) - (c/(4*R))* ( (lambda + (1-a).*sin(phi)).^2 +...
         (1-a).^2.* cos(phi).^2).*cla.*(atan((1-a).*cos(phi))./(lambda + ...
         (1-a).*sin(phi))).*(cos(phi+(atan((1-a).*cos(phi))./(lambda +...
         (1-a).*sin(phi))))./cos(phi)); %equazione in a di cui trovare lo zero
    a(1,i) = (fzero(eq,0.05)); % zero della funzione, per ogni phi abbiamo un valore di a --> a è costante lungo l'ala ma varia con phi
    ind(j,:) = a; %lungo la riga varia con phi, lungo la colonna varia con lambda
end

v_vinf(j,:) = sqrt((lambda + (1-a).*sin(phiv)).^2 + (1-a).^2 .*cos(phiv).^2); %v/vinf: lungo le righe varia con phi, lungo le colonne varia con lambda
alpha(j,:) = (atan((1-a).*cos(phiv))./(lambda + (1-a).*sin(phiv))); %alpha: lungo le righe varia con phi, lungo le colonne varia con lambda
cost_p = (N*c*lambda)./(4*pi*R); % costante nel calcolo dell'integrale (cp)
cp(j,1) = cost_p.* trapz(phiv, (v_vinf(j,:).^2*cla.*alpha(j,:).*sin(alpha(j,:)).*(1-(cd./(cla.*alpha(j,:))).*cot(alpha(j,:))))); %calcolo di cp
% cost_q = (N*c)./(8*pi^2);
% cq(j,1) = cost_q.* trapz(phiv, (v_vinf*cla.*alpha.*sin(alpha).*(1-(cd./(cla.*alpha)).*cot(alpha))));
cq(j,1) = cp(j,1)/lambda; %9.6 pagina 118 %calcolo cq
%cp(j,1).*R./lambda;
end

%% Plot
figure(1)
plot(lambdav,cp','k.-'); grid on;
xlabel('\lambda'); ylabel('C_P');
%cp è per unità di lunghezza

figure(2)
plot(lambdav,cq','k.-'); grid on;
xlabel('\lambda'); ylabel('C_Q');
%cq per unità di lunghezza

figure(3)
plot(lambdav,ind(:,1),'k.-') 
xlabel('\lambda'); ylabel('a')
legend('(\phi = 0)','Location','best')
title('induzione a al variare di lambda')
grid on

figure(4)
plot(phiv,ind(1,:),'k.-')
xlabel('\phi'); ylabel('a');
legend('(\lambda = 2)','Location','best')
title('induzione a al variare di phi')
grid on


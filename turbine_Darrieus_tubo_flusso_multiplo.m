function [cq,cp] = turbine_Darrieus_tubo_flusso_multiplo (alphamax,c,R,cla,N,cd)
phiv = linspace(0,360,100); %dominio di phi (posizione angolare dell'ala)
phiv = deg2rad(phiv); % dominio di phi in rad
lambdav = linspace(2,17.5,100); %vettore dei lambda
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
          a(1,i) = (fzero(eq,0.2)); % zero della funzione, per ogni phi abbiamo un valore di a --> a è costante lungo l'ala ma varia con phi
          ind(j,:) = a; %lungo la riga varia con phi, lungo la colonna varia con lambda

    v_vinf(j,:) = sqrt((lambda + (1-a).*sin(phiv)).^2 + (1-a).^2 .*cos(phiv).^2); %v/vinf: lungo le righe varia con phi, lungo le colonne varia con lambda
    alpha(j,:) = atan2(((1-a).*cos(phiv)),(lambda + (1-a).*sin(phiv))); %alpha: lungo le righe varia con phi, lungo le colonne varia con lambda
    cost_p = (N*c*lambda)./(4*pi*R); % costante nel calcolo dell'integrale (cp)
    cp(j,1) = cost_p.* trapz(phiv, (v_vinf(j,:).^2*cla.*alpha(j,:).*sin(alpha(j,:)).*(1-(cd./(cla.*alpha(j,:))).*cot(alpha(j,:))))); %calcolo di cp
    cq(j,1) = cp(j,1)/lambda; %9.6 pagina 118 %calcolo cq
    end

massimi = zeros(numel(lambdav),1);
alphadeg = rad2deg(alpha); % dominio di phi in rad

    for h = 1 : numel(phiv)
        massimi(h) = max(alphadeg(h,:));
    end

contatorelambdamin = 0;

    for f = 1 : numel(phiv)
        if alphamax < massimi(f)
            contatorelambdamin = contatorelambdamin + 1;
        end
    end

contatorelambdamax = 0;

    for k = 1 : length(phiv)
        if cp(k) < 0 
        contatorelambdamax = contatorelambdamax + 1;
        end
    end

end
% cq e cp vengono ridimensionati considerando il
% lambdaminimo(alpha deve essere minore dell'alphadistallo) e il lambdamassimo(cp>0)
cq = cq(contatorelambdamin:end-contatorelambdamax,1);
cp = cp(contatorelambdamin:end-contatorelambdamax,1);
end




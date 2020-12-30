%% test case with Weick 1931 (pag 166 figura C-2)
clear all; close all; clc
V= 40; %m/s
rpm= 4500; %giri al minuto
R_dim= 1.58; %m
D=1.58*2; %m
n=rpm/60; %giri al secondo
%R= 1
r_hub_dim= 0.14; %m
r_hub=r_hub_dim/R_dim;
N= 3; %numero di pale
ns= 10; %numero sezioni
r= linspace(r_hub_dim, 0.97*R_dim, ns);
r= r/1.58;
R=R_dim/1.58;
c= 0.01*[9, 15, 17, 18.5, 19, 19, 17.5, 15.5, 12.5, 9]; %m
Cl_alpha= ones(1, ns);
Cl_alpha= 6.28*Cl_alpha;
beta_deg= [40, 39.7, 37.6, 32.3, 27.8, 24.4, 21.8, 19.9, 18.3, 17.6];
x=linspace(r_hub,R,100);

%omega= convangvel(rpm, 'rpm', 'rad/s')
omega= rpm*2*pi/60;
i=1;
while i< ns+1
beta(i)= convang(beta_deg(i), 'deg', 'rad');
i= i+1;
end

%% Metodo A

lam= V/(omega*R);
sigma= N*c./(pi*r);
Vt= omega*r;
Vr= Vt.*sqrt((lam^2)+ (r.^2));
phi= atan ((V)./(Vt));

alpha_i_A= 0.5*(-(lam./r + sigma.*Cl_alpha./(8*r.^2).*(Vr./Vt)) +...
    sqrt((lam./r + sigma.*Cl_alpha./(8*r.^2).* (Vr./Vt)).^2 +...
    sigma.*Cl_alpha.*Vr.*(beta-phi)./(2*(r.^2).*Vt)));
alpha_i_A_spline=interp1(r,alpha_i_A,x,'spline');

%% Metodo B
lam= V/(omega*R);
sigma= N*c./(pi*r);
Vt= omega.*r;
Vr= Vt.*sqrt((lam^2)+ (r.^2));
phi= atan ((V)./(Vt));

F=(2/pi)*acos(exp(N/(2*lam)*(r-R)/R)); %fattore correttivo di Prandtl

alpha_i_B= 0.5*(-(lam./r + sigma.*Cl_alpha./(8*r.*F.*cos(phi))) ...
         + sqrt((lam./r + sigma.*Cl_alpha./(8*r.*F.*cos(phi))).^2 ...
         + sigma.*Cl_alpha.*(beta-phi)./(2*r.*F.*cos(phi))));
alpha_i_B_spline=interp1(r,alpha_i_B,x,'spline');

w_t_B= Vr.*alpha_i_B.*sin(phi+alpha_i_B); % formula C-8
w_t_C= sigma.*Cl_alpha./(8*r.*F).*Vr.*(beta- phi-alpha_i_B); %%
w_a_B= Vr.*alpha_i_B.*cos(phi+alpha_i_B); % formula C-9

%% Metodo C (con Newton-Raphson) 

alfaisa= alpha_i_B;
wtsa= w_t_B;
wasa= w_a_B;

w_t_B= w_t_B-(0.1*w_t_B);

%faccio partire la prima iterazione prima di entrare nel ciclo

test1= V^2+ (4*w_t_B.*((omega*R*r)- w_t_B)); %radice della C-15
test2= lam^2+ (4*w_t_B./(Vt)).*(r-(w_t_B./(Vt))); %radice della C-17

 if (test1 <0 | test2 <0)
     w_t_B=wtsa;
     w_a_B=wasa;
 else
     
term1= (sqrt((V^2)+(4*w_t_B.*(omega*R*r- w_t_B)))-V); % 2w_a_C dalla C-15
term2= lam+ sqrt((lam^2)+(4*(w_t_B./Vt).*(r-(w_t_B./Vt)))); %seconda radice di C-17

aa= beta-(atan((w_t_B*2)./(term1))); % termine in parentesi nella C-14
bb= sqrt((0.25*(term2.^2)+((r-w_t_B./Vt)).^2)); %tutto il termine sotto radice della C-17
cc= 8*r.*F.*w_t_B./Vt; % ha portato al primo membro tutti i termini della C-17 tranne sigma,
                      % Cl_alpha, aa, bb

da= ((-2)./(((term1).^2)+ 4*w_t_B.*w_t_B)).*(term1...
    -(w_t_B.*(((2*omega*R*r)-(4*w_t_B))./(term1+V))));      %derivata di aa

db= (((term2.*((r./Vt)-((2*w_t_B)./((Vt).^2))))...
    ./(2*(term2- lam)))+ (w_t_B./((Vt).^2))-(r./Vt))./bb;   %derivata di bb

dc= 8*r.*F./Vt;  %derivata di cc

y= (sigma.*Cl_alpha.*aa.*bb)- cc; % y indica la differenza tra cc e..
                                  % il resto del termine della C-17
dy= (sigma.*Cl_alpha.*((aa.*db)+(bb.*da)))-dc;

i=1;

epsi= 1.e-2; % nel listato fortran vi è un termine 'err' che ...
             % non è specificato viene preso dall'esterno 
if abs(y./dy) < epsi
w_t_B= w_t_B- (y./dy);
else
    while ((abs(y./dy)> epsi) | i==1001)
    
    test1= V^2+ (4*w_t_B.*((omega*R*r)- w_t_B));
    test2= lam^2+ (4*w_t_B./(Vt)).*(r-(w_t_B./(Vt)));

 if (test1 <0 | test2 <0)
     w_t_B=wtsa;
     w_a_B=wasa;
 else
        
term1= (sqrt((V^2)+(4*w_t_B.*(omega*R*r- w_t_B)))-V); % 2w_a_C dalla C-15
term2= lam+ sqrt((lam^2)+(4*(w_t_B./Vt).*(r-(w_t_B./Vt)))); %seconda radice di C-17

aa= beta-(atan((w_t_B*2)./(term1))); % termine in parentesi nella C-14
bb= sqrt((0.25*(term2.^2)+((r-w_t_B./Vt)).^2)); %tutto il termine sotto radice della C-17
cc= 8*r.*F.*w_t_B./Vt; % ha portato al primo membro tutti i termini della C-17 tranne sigma,..
                      % Cl_alpha, aa, bb

da= ((-2)./(((term1).^2)+ 4*w_t_B.*w_t_B)).*(term1...
    -(w_t_B.*(((2*omega*R*r)-(4*w_t_B))./(term1+V))));

db= (((term2.*((r./Vt)-((2*w_t_B)./((Vt).^2))))...
    ./(2*(term2- lam)))+ (w_t_B./((Vt).^2))-(r./Vt))./bb;

dc= 8*r.*F./Vt;

y= (sigma.*Cl_alpha.*aa.*bb)- cc;

dy= (sigma.*Cl_alpha.*((aa.*db)+(bb.*da)))-dc;

i= i+1;

w_t_B= w_t_B-(0.1*w_t_B);

 end
  
    end
    q= y./dy;
        w_t_B= w_t_B- (y./dy);
        w_a_C= 0.5*(sqrt(V.^2+ 4*w_t_B.*(omega*R*r- w_t_B))- V);
        alpha_i_C= atan (w_t_B./ w_a_C)- phi;
        alpha_i_C_spline=interp1(r,alpha_i_C,x,'spline');
end
        
   if i>1000
        disp('non ci sono sol per wt')
   end
end
 
%% Calcolo  dei coefficienti aerodinamici
J=V/(n*D);
Cd=0.01;

delta_ct_a=(pi/8)*sigma.*(J^2 + (pi*r).^2).*(Cl_alpha.*alpha_i_A.*(cos(phi+alpha_i_A))-...
            Cd*sin(phi+alpha_i_A));
delta_ct_b=(pi/8)*sigma.*(J^2 + (pi*r).^2).*(Cl_alpha.*alpha_i_B.*(cos(phi+alpha_i_B))-...
            Cd*sin(phi+alpha_i_B));
delta_ct_c=(pi/8)*sigma.*(J^2 + (pi*r).^2).*(Cl_alpha.*alpha_i_C.*(cos(phi+alpha_i_C))-...
            Cd*sin(phi+alpha_i_C));
        
delta_cp_a=((r*pi^2)/8).*sigma.*(J^2 + (pi*r).^2).*(Cl_alpha.*alpha_i_A.*(sin(phi+alpha_i_A))+...
            Cd*cos(phi+alpha_i_A));
delta_cp_b=((r*pi^2)/8).*sigma.*(J^2 + (pi*r).^2).*(Cl_alpha.*alpha_i_B.*(sin(phi+alpha_i_B))+...
            Cd*cos(phi+alpha_i_B));
delta_cp_c=((r*pi^2)/8).*sigma.*(J^2 + (pi*r).^2).*(Cl_alpha.*alpha_i_C.*(sin(phi+alpha_i_C))+...
            Cd*cos(phi+alpha_i_C));

delta_ct_a_spline=interp1(r,delta_ct_a,x,'spline');
delta_ct_b_spline=interp1(r,delta_ct_b,x,'spline');
delta_ct_c_spline=interp1(r,delta_ct_c,x,'spline');
delta_cp_a_spline=interp1(r,delta_cp_a,x,'spline');
delta_cp_b_spline=interp1(r,delta_cp_b,x,'spline');
delta_cp_c_spline=interp1(r,delta_cp_c,x,'spline');

%% Grafici 
figure
plot(x, alpha_i_A_spline, 'k')
xlabel('raggio adim');ylabel('\alpha_i [rad]');
title ('metodo A')

figure
plot(x, alpha_i_A_spline, 'k')
hold on
plot(x, alpha_i_B_spline, 'b')
xlabel('raggio adim');ylabel('\alpha_i [rad]');
legend('Metodo A','Metodo B')
title ('metodo A e B a confronto')

figure
plot(x, alpha_i_A_spline, 'k')
hold on
plot(x, alpha_i_B_spline, 'b')
hold on
plot(x, alpha_i_C_spline, 'r')
xlabel('raggio adim');ylabel('\alpha_i [rad]');
legend('Metodo A','Metodo B','Metodo C')
title ('i 3 metodi a confronto')

figure
plot(x,delta_ct_a_spline,'k',x,delta_ct_b_spline,'b',x,delta_ct_c_spline,'r')
xlabel('raggio adim');ylabel('dC_T/dr');
legend('Metodo A','Metodo B','Metodo C')

figure
plot(x,delta_cp_a_spline,'k',x,delta_cp_b_spline,'b',x,delta_cp_c_spline,'r')
xlabel('raggio adim');ylabel('dC_P/dr');
legend('Metodo A','Metodo B','Metodo C')
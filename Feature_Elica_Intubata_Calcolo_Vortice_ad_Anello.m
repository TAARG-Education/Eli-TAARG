function [f] = Feature_Elica_Intubata_Calcolo_Vortice_ad_Anello(ch,Di,De)
%Funzione per il calcolo della velocità indotta da un vortice ad anello. Si
%utilizza la legge di Biot-Savart per il calcolo delle tre componenti,
%scritta in un riferimento cilindrico.
%La legge di Biot-Savart così scritta è: data da
%                 V(P)=-(Gamma/4*pi)*integ(0,2pi)((R x dl)/R^3).
%Le tre componenti dimensionali della legge di Biot-Savart nelle tre
%direzioni vengono così definite:
%       V_x=(Gamma/4*pi)*r*z*integ(0,2pi)((cos(theta)dtheta)/R^3)   dir x
%       V_y=-(Gamma/4*pi)*r*z*integ(0,2pi)((sin(theta)dtheta)/R^3)  dir y
%V_z=-(Gamma/4*pi)*r*integ(0,2pi)(((x-x^)cos(theta)+(y-y^)sin(theta)dtheta)/R^3)

%Si adimensionalizzano le formule per una velocità caratteristica che è la
%circolazione sul diametro del vortice ad anello definita come
%Gamma*Dvortice.
%Calcoliamo di base le tre quantità integrali

%Gui per l'input da parte dell'utente.
%  title='Dati Input';
% prompt={'Lunghezza carenatura[m]','Diametro ingresso [m]','Diametro uscita[m]'};
% answer=inputdlg(prompt,title,[1 55;1 55;1 55]);
% ch=str2double(answer{1});
% Di=str2double(answer{2});
% De=str2double(answer{3});

ch=3; %m  
Di=4.5; %m
De=1.20; %m

% [f1,c_D14,RatioD,D14,D34] = Universale2(Di,De,ch);
t=linspace(0,1,100);             %commentare da riga 31 a 40 se si utilizza la funzione 30 per il calcolo dei rapporti c_D14 e RatioD 
for i=1:length(t) 
xs(i)=0+t(i)*(ch);
ys(i)=Di/2+t(i)*(De/2-Di/2);
end
D14=2*ys(25);
D34=2*ys(75);

c_D14=ch/D14 ;
RatioD=D34/D14;

r=(D14/2);
x=D34/2;
z=ch*(3/4);
y=0;

R1=@(t) cos(t)./(((x-r.*cos(t)).^2 +((y-r.*sin(t)).^2)+z^2)).^(3/2);
R2=@(t) sin(t)./(((x-r.*cos(t)).^2 +(y-r.*sin(t)).^2+z^2)).^(3/2);
R3=@(t) ((x-r.*cos(t)).*cos(t)+(y-r.*sin(t)).*sin(t))./(((x-r.*cos(t)).^2 +(y-r.*sin(t)).^2+z^2)).^(3/2);
Ix=integral(R1,0,2*pi);
Iy=integral(R2,0,2*pi);
Iz=integral(R3,0,2*pi);
%Adimensionalizzazione;
fx=(r*z*D34/4)*Ix;
fy=(r*z*D34/4)*Iy;
fz=(r*D34/4)*Iz;
%Somma vettoriale;
f=sqrt(fx^2+fy^2+fz^2);

% P=abs(f1/f-1)*100;
%Grafico in cui ottenere i dati:
figure(2);
plot(f,f)
axis ([0 1 0 1]);
text(0,1,'Velocità Indotta');
% text(0,0.90,['Metodo grafico, f= ',num2str(f1)]);
text(0,0.80,['Metodo matematico, f= ',num2str(f)]);
% text(0,0.70,['Errore percentuale, %= ',num2str(P)]);
axis off;
end


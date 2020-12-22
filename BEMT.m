%% BEMT_Function Summary
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
% |Name        : Function BEMT (A, B and C Methods)                                           |
% |Author      : Crispino Umberto & Di Maio Daniele                                           |
% |              University of Naples Federico II.                                            |
% |Version     : 1.0.0                                                                        |
% |Date        : 15/12/2020                                                                   |
% |Modified    : 22/12/2020                                                                   |
% |Description : Using flight conditions and the geometry of the propeller, through three     | 
%                different methods the induced angle of incidence is obtained                 |
% |Reference   :                                                                              |
% |Input       : Geometry of the propeller (number of blades (N), radius (R), number of       |
%                sections (ns) and radius of the hub) and flight condintions (V and rpm)      |  
% |Output      : Induced incidence angle valued by A, B and C Methods                         |
% |Note        :                                                                              |
% ==============================================================================================

function [alpha_i_A, alpha_i_B, alpha_i_C, beta, c, Cl_alpha]=BEMT(V, rpm, R, r_hub, N, ns)

i=1;
while i< ns+1
s=input(['inserisci valore della corda (m) alla stazione  ',num2str(i), ' : ']);
c(i)=s;
i = i+1;
end
i=1;
while i< ns+1
f=input(['inserisci valore del Cl_alpha (1/rad) alla stazione  ',num2str(i), ' : ']);
Cl_alpha(i)=f;
i = i+1;
end
i=1;
while i< ns+1
d=input(['inserisci valore dell angolo di calettamento aerodinamico (deg) alla stazione  ',num2str(i), ' : ']);
beta_deg(i)=d;
i=i+1;
end

omega=convangvel(rpm/60, 'rpm', 'rad/s');
i=1;
while i< ns+1
beta(i)=convang(beta_deg(i), 'deg', 'rad');
i= i+1;
end
lam=V/(omega*R);
r=linspace(r_hub, R, ns);
sigma=N*c./(2*pi*r); %solidità della stazione della pala
Vt=omega.*r;
Vr=Vt.*sqrt((lam^2)+(r.^2));
phi=atan2((V),(Vt));  
F=(2/pi)*acos(exp(N/(2*lam)*(r-R)/R)); %fattore correttivo di Prandtl

%% A-Method
alpha_i_A= -0.5*(-(lam./r + sigma.*Cl_alpha./(8*r.^2).*(Vr./Vt)) ...
         + sqrt((lam./r + sigma.*Cl_alpha./(8*r.^2).*(Vr./Vt)).^2 ...
         + sigma.*Cl_alpha.*Vr.*(beta-phi)./(2*(r.^2).*Vt)));

%% B-Method
alpha_i_B= 0.5*(-(lam./r + sigma.*Cl_alpha./(8*r.*F.*cos(phi))) ...
         + sqrt((lam./r + sigma.*Cl_alpha./(8*r.*F.*cos(phi))).^2 ...
         + sigma.*Cl_alpha.*(beta-phi)./(2*r.*F.*cos(phi))));

%% C-Method
w_t_B= sigma.*Cl_alpha./(8*r.*F).*Vr.*(beta- phi-alpha_i_B);

w_a_C= 0.5*(sqrt(V.^2+ 4*w_t_B.*(omega*R*r- w_t_B))- V);

Ve2= (omega*R*r- w_t_B).^2+ (V+ 0.5*sqrt(V.^2+ 4*w_t_B.*(omega*R*r- w_t_B)-0.5*V)).^2;
Ve=sqrt(Ve2);

w_t_C= ((sigma.*Ve.*Cl_alpha)./(8*r.*F)).*(beta-atan2(w_t_B,w_a_C))...
       .*sqrt((r- w_t_B/(omega*R)).^2+ ...
       0.25*(sqrt(lam.^2+ 4*w_t_B.*(r- w_t_B/(omega*R))/(omega*R))+ lam).^2);

w_t_C_old= w_t_C;
w_t_C_new= ((sigma.*Ve.*Cl_alpha)./(8*r.*F)).*(beta-atan2(w_t_C_old, w_a_C))...
       .*sqrt((r- w_t_C_old/(omega*R)).^2+ ...
       0.25*(sqrt(lam.^2+ 4*w_t_C_old.*(r- w_t_C_old/(omega*R))/(omega*R))+ lam).^2);

w_a_C= 0.5*(sqrt(V.^2+ 4*w_t_C_new.*(omega*R*r- w_t_C_new))- V);

Ve2= (omega*R*r- w_t_C_new).^2+ ...
     (V+ 0.5*sqrt(V.^2+ 4*w_t_C_new.*(omega*R*r- w_t_C_new)-0.5*V)).^2;
Ve=sqrt(Ve2);
   
   
epsi= 1.e-3;

while abs(w_t_C_new-w_t_C_old) > epsi

    w_t_C_old= w_t_C_new; %old i= new i-1

    w_a_C= 0.5*(sqrt(V.^2+ 4*w_t_C_old.*(omega*R*r- w_t_C_old))- V);

    Ve2= (omega*R*r- w_t_C_old).^2+ ...
         (V+ 0.5*sqrt(V.^2+ 4*w_t_C_old.*(omega*R*r- w_t_C_old)-0.5*V)).^2;
    Ve=sqrt(Ve2);

    w_t_C_new= ((sigma.*Ve.*Cl_alpha)./(8*r.*F)).*(beta-atan2(w_t_C_old, w_a_C))...
       .*sqrt((r- w_t_C_old/(omega*R)).^2+ ...
       0.25*(sqrt(lam.^2+ 4*w_t_C_old.*(r- w_t_C_old/(omega*R))/(omega*R))+lam).^2); %new i

end
alpha_i_C= atan2 (w_t_C_new, w_a_C) - phi;

end


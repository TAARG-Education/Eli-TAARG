%% BEMT_function Summary
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
% |Name        : Function BEMT (1, 2 and 3 Methods)                                           |
% |Author      : Crispino Umberto & Di Maio Daniele                                           |
% |              University of Naples Federico II                                             |
% |Version     : 1.0.0                                                                        |
% |Date        : 15/12/2020                                                                   |
% |Modified    : 26/02/2021                                                                   |
% |Description : Using flight conditions and the geometry of the propeller, through three     | 
%                different methods, the induced angle of incidence is obtained                |
% |Reference   : "Lezioni di Aerodinamica dell'ala rotante" by Renato Tognaccini              |
% |Input       : Geometry of the propeller (number of blades (N), radius (R), number of       |
%                sections (ns) and radius of the hub (r_hub)) and flight conditions           |
%                (V and rpm)                                                                  |  
% |Output      : Induced incidence angle valued by Methods N1, N2 and N3                      |
% |Note        :                                                                              |
% ==============================================================================================

function [alpha_i_deg, beta_deg, c, Cl_alpha]=BEMT(V, rpm, R, r_hub, N, ns)

p = 1;
while p < ns+1
s=input(['INSERT THE VALUE OF THE CHORD [m] AT THE STATION ',num2str(p), ' : ']);
c(p)=s;
p = p+1;
end
p = 1;
while p < ns+1
f=input(['INSERT THE VALUE OF THE Cl_ALPHA [1/rad] AT THE STATION ',num2str(p), ' : ']);
Cl_alpha(p)=f;
p = p+1;
end
p = 1;
while p < ns+1
d=input(['INSERT THE VALUE OF THE PITCH ANGLE [deg] AT THE STATION ',num2str(p), ' : ']);
beta_deg(p)=d;
p = p+1;
end

omega  = convangvel(rpm, 'rpm', 'rad/s');
r      = linspace(r_hub, R*0.99, ns);
beta   = convang(beta_deg, 'deg', 'rad');
r_ndim = r/R;
x      = linspace(r_hub,R,100);

lam    = V/(omega*R);
sigma  = N*c./(pi*R);
Vt     = omega*R;
Vr     = sqrt((V^2)+ (omega*R.*r_ndim).^2);
phi    = atan((V)./(omega*r)); 

disp('Which theory do you want to use?')
disp('1 - IMPULSIVE THEORY')
disp('2 - WHIRLING THEORY WITH SMALL DISTURBANCES')
disp('3 - WHIRLING THEORY')

l= input(['']);
    
%% METHOD N1
if l==1
    
alpha_i = 0.5*(-(lam./r_ndim + sigma.*Cl_alpha.*(Vr./Vt)./(8*r_ndim.^2)) +...
          sqrt((lam./r_ndim + sigma.*Cl_alpha.*(Vr./Vt)./(8*r_ndim.^2)).^2 +...
          sigma.*Cl_alpha.*Vr.*(beta-phi)./(2*(r_ndim.^2).*Vt)));

alpha_i_spline = interp1(r,alpha_i,x,'spline');
else
%% METHOD N2
if l==2
    
   if phi(ns)== 0
   F= ones (1, ns);
   else       
   F=(2/pi)*acos(exp(N/(2*lam)*((r-R)/R)));
   end
   
alpha_i = 0.5*(-(lam./r_ndim + sigma.*Cl_alpha./(8*r_ndim.*F.*cos(phi))) + ...
          sqrt((lam./r_ndim + sigma.*Cl_alpha./(8*r_ndim.*F.*cos(phi))).^2 + ...
          sigma.*Cl_alpha.*(beta-phi)./(2*r_ndim.*F.*cos(phi))));

alpha_i_spline = interp1(r,alpha_i,x,'spline');
else
%% METODO 3
if l==3
   max_iter=1000;
   warndlg('WARNING: questo metodo è da rivedere ai posteri... o forse no');
    
       if phi(ns)== 0
       F= ones (1, ns);
       else       
       F=(2/pi)*acos(exp(N/(2*lam)*((r-R)/R)));
       end
    
alpha_i_B = 0.5*(-(lam./r_ndim + sigma.*Cl_alpha./(8*r_ndim.*F.*cos(phi))) ...
          + sqrt((lam./r_ndim + sigma.*Cl_alpha./(8*r_ndim.*F.*cos(phi))).^2 ...
          + sigma.*Cl_alpha.*(beta-phi)./(2*r_ndim.*F.*cos(phi))));

w_t_B = Vr.*alpha_i_B.*sin(phi+alpha_i_B); % formula C-8
w_t_C = sigma.*Cl_alpha./(8*r.*F).*Vr.*(beta- phi-alpha_i_B); 
w_a_B = Vr.*alpha_i_B.*cos(phi+alpha_i_B); % formula C-9

alfaisa = alpha_i_B;
wtsa    = w_t_B;
wasa    = w_a_B;
w_t_B   = w_t_B-(0.1*w_t_B);
value   = 10* ones(1, ns);      
epsi    = 1e-10; 
        
for j= 1:ns
    i=0;
if value(j) < epsi
    w_t_B(j)= w_t_B(j)-(y(j)/dy(j));
else
    
while [value(j)> epsi , i<max_iter+1]
    i= i+1;
    test1(j)= V^2+ (4*w_t_B(j).*((omega*R*r_ndim(j))- w_t_B(j))); %radice della C-15         
    test2(j)= lam^2+ (4*w_t_B(j)./(Vt)).*(r_ndim(j)-(w_t_B(j)./(Vt))); %radice della C-17    

if (test1 <0 | test2 <0)
     w_t_B(j)=wtsa(j);
     w_a_B(j)=wasa(j);
else

term1(j)= (sqrt(test1(j))-V);  %2*w_a_C dalla C-15          
term2(j)= lam+ sqrt(test2(j)); %seconda radice di C-17       

aa(j)= beta(j)-(atan((w_t_B(j)*2)./(term1(j)))); % termine in parentesi nella C-14
bb(j)= sqrt((0.25*(term2(j).^2)+((r_ndim(j)-w_t_B(j)./Vt)).^2)); %tutto il termine sotto radice della C-17
cc(j)= 8*r_ndim(j).*F(j).*w_t_B(j)./Vt; % ha portato al primo membro tutti i termini della C-17 tranne sigma,..
                      % Cl_alpha, aa, bb 

da(j)= ((-2)./(((term1(j)).^2)+ 4*w_t_B(j).*w_t_B(j))).*(term1(j)... %derivata di aa
     -(w_t_B(j).*(((2*omega*R*r_ndim(j))-(4*w_t_B(j)))./(term1(j)+V))));
db(j)= (((term2(j).*((r_ndim(j)./Vt)-((2*w_t_B(j))./((Vt).^2))))...  %derivata di bb
    ./(2*(term2(j)- lam)))+ (w_t_B(j)./((Vt).^2))-(r_ndim(j)./Vt))./bb(j);
dc(j)= 8*r_ndim(j).*F(j)./Vt;  %derivata di cc

y(j)= (sigma(j).*Cl_alpha(j).*aa(j).*bb(j))- cc(j);
dy(j)= (sigma(j).*Cl_alpha(j).*((aa(j).*db(j))+(bb(j).*da(j))))-dc(j);

value(j)= abs(y(j)./dy(j));
w_t_B(j)= w_t_B(j)-(y(j)./dy(j));

alpha_i(j)= beta(j)-aa(j)-phi(j);

end
end
end
I(j)=i;
end

alpha_i_spline=interp1(r,alpha_i,x,'spline');
    
if i>max_iter
   disp('There is no value for w_t');
end
end
end
end

alpha_i_deg        = convang(alpha_i,'rad','deg');
alpha_i_deg_spline = convang(alpha_i_spline,'rad','deg');

%% Chart
figure
plot(x, alpha_i_deg_spline, 'k')
xlabel('r [m]'); ylabel('\alpha_i [deg]');
axis ([0 R min(alpha_i_deg)-1 max(alpha_i_deg)+1]);
end
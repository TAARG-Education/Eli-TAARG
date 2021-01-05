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
% |Name        : turbine_Darrieus_tubo_flusso_singolo.m                                       |
% |Author      : Marco Gugliucciello - Angelo Dasco                                    |
% |              University of Naples Federico II.                                             |
% |Version     : 1.5                                                                           |
% |Date        : 05/06/2020                                                                    |
% |Modified    :                                                                     |
% |Description : characteristic curves evaluation for a Darrieus turbine through the           |
%                 multiple streamtubes theory                                                  |
% |Reference   : Tognaccini R., (2019), "Lezioni di Aerodinamica dell'ala rotante"             |
%                De Vries O., (1979), "Fluid Dynamic Aspects of Wind Energy Conversion"        |
% |Input       :  (alphamax) = angle of attack at the stall of the considered airfoil          |
%                 (c) = chord of the blade                                                     |
%                 (R) = blade radius                                                           |
%                 (Cla) = slope coefficient of the Cl curve for the considered airfoil         |
%                 (N) = number of blades                                                       |
%                 (Cd) = drag coefficient of the considered airfoil                            |
% |Output      :  (cp) = vector of power coefficient values for the                            |
%                 turbine                                                                      |
%                 (cq) = vector of torque coefficient values for the                           |
%                 turbine                                                                      |
%                 (lambdav) = vector of tip speed values corresponding to                      |
%                 the cp and cq coefficients                                                   |
% |Note        :                                                                               |
% ==============================================================================================

%Create a domain vector. 
        Cla=6.28;
        alphamax = 14;
        c = 1;
        R = 30;
        cla = 6.28;
        N = 3;
        cd = 0.01;
           
y1=linspace(0.01,20,100);
phiv = linspace(0,360,100);
phiv = deg2rad(phiv);
ind = zeros(1,length(phiv));

%equations


   
for i = 1:length(y1)
    
    y=y1(i);
    
    Y= ((y-sin(phiv)).^2+cos(phiv).^2).*Cla.*atan(cos(phiv)./(y+sin(phiv))).*cos(phiv+atan(cos(phiv)./(y+sin(phiv))));
    


% Integrate Y using trapz. 
A =   trapz(phiv,Y);
x = 1 + (N.*c./(8.*R.*2.*pi)).*A;
a(i)= 1-1/x;
lambda(i) = y.*(1-a(i));



end
v_vinf = sqrt((lambda + (1-a).*sin(phiv)).^2 + (1-a).^2 .*cos(phiv).^2);
        alpha = atan2(((1-a).*cos(phiv)),(lambda + (1-a).*sin(phiv)));
        cost_p = (N*c*lambda)./(4*pi*R);
        % computing Cp - numerical integration
        cp = cost_p.* trapz(phiv, (v_vinf.^2*cla.*alpha.*sin(alpha).*...
            (1-(cd./(cla.*alpha)).*cot(alpha))));
        cq = cp/lambda; % Cq value

        plot(lambda,cp,'k-');
        legend('single streamtube theory','Location','best');
        h=xlabel('$\lambda$','Interpreter','latex'); set(h,'Fontsize',20);
        h=ylabel('$C_P$','Interpreter','latex');
        set(h,'Fontsize',20);
        title('\sigma = 0.1 - Cd = 0.01');

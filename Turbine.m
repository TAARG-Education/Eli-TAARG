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
% |Name        : turbine_Darrieus_tubo_flusso_singolo.m                                        |
% |Author      : Marco Gugliucciello - Angelo Dasco                                            |
% |              University of Naples Federico II.                                             |
% |Version     : 2                                                                             |
% |Date        : 05/06/2020                                                                    |
% |Modified    : 15/01/2021                                                                    |
% |Description : characteristic curves evaluation for a Darrieus turbine through the           |
%                 multiple streamtubes theory                                                  |
% |Reference   : Tognaccini R., (2019), "Lezioni di Aerodinamica dell'ala rotante"             |
%                De Vries O., (1979), "Fluid Dynamic Aspects of Wind Energy Conversion"        |
% |Input       :  (c) = chord of the blade                                                     |
%                 (R) = blade radius                                                           |
%                 (Cla) = slope coefficient of the Cl curve for the considered airfoil         |
%                 (N) = number of blades                                                       |
%                 (Cd) = drag coefficient of the considered airfoil                            |
%                                                                                              |
% |Output      :  (cp) = vector of power coefficient values for the                            |
%                 turbine                                                                      |
%                 (cq) = vector of torque coefficient values for the                           |
%                 turbine                                                                      |
%                 (lambda) = vector of tip speed values corresponding to                       |
%                 the cp and cq coefficients                                                   |
% |Note        :  The function works                                                           |
% ==============================================================================================
function [cp2, lambda2] = Darrieus_singolo(Cla, c, R, N, cd,alphamax)

%Domain vector, initial values. 
y1=linspace(0.01,20,100);       %Introduction of the variable y
phiv = linspace(0,360,100);     %Phi domain (deg)
phiv = deg2rad(phiv);           %Phi domain (rad)  

%Equations 
 
for i = 1:length(y1)            % a, lambda, cp, cq computing 
    
    y=y1(i);
    
%   Computing of the induction factor a

    Y= ((y-sin(phiv)).^2+cos(phiv).^2).*Cla.*atan(cos(phiv)./(y+sin(phiv))).*...
        cos(phiv+atan(cos(phiv)./(y+sin(phiv))));   
    A =   trapz(phiv,Y);
    x = 1 + (N.*c./(8.*R.*2.*pi)).*A;
    a= 1-1/x; 

%   Determination of tip speed ratio

    lambda = y.*(1-a);
    lambda2(i,1)=lambda;

%   Determination of velocity ratio
      v_vinf = sqrt((lambda + (1-a).*sin(phiv)).^2 +...
          (1-a).^2 .*cos(phiv).^2);

%   Determination of angle of attack 
       alpha = atan2(((1-a).*cos(phiv)),(lambda + (1-a).*sin(phiv)));    
      
%   Computing Cp - numerical integration
       cost_p = (N*c*lambda)./(4*pi*R);            % Costant part 
       cpint =trapz(phiv,(v_vinf.^2.*Cla.*alpha.*sin(alpha).*...
           (1-(cd./(Cla.*alpha)).*cot(alpha))));   % Integral
        
%       
       cp=cost_p.*cpint;                           % Cp value
       cq = cp/lambda;                             % Cq value
        
       cp2(i,1)=cp;
end

    max_v = zeros(numel(lambda),1);
    alphadeg = rad2deg(alpha); % [deg]

%   Find max alpha for each row of the matrix alpha
for h = 1 : numel(phiv)
    max_v(h) = max(alphadeg(h));
end

    contatorelambdamin = 0; %intial value for lambda min index

%c  heck on stall angle
for f = 1 : numel(phiv)
    if alphamax < max_v(f)
        contatorelambdamin = contatorelambdamin + 1;
    end
end

end



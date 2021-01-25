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
% |Name        : Curve_caratteristiche_aeromotori_asse_orizzontal.m                            |
% |Author      : Nasti Giuseppe, Tomasso Armando Diego                                         |
% |              University of Naples Federico II.                                             |
% |Version     : 1.0                                                                             |
% |Date        : 20/01/2021                                                                    |
% |Modified    : 20/01/2021                                                                    |
% |Description : characteristic curves evaluation for a windmill                               |
%                                                                                              |
% |Reference   : Tognaccini R., (2019), "Lezioni di Aerodinamica dell'ala rotante"             |
%                BURTON, T., JENKINS, N.,SHARPE, D., BOSSANYI, E.,(2011),                      |
%                             "Wind Energy Handbook", 2nd Edition, Wiley and Sons.             |
% |Input       :  (c) = chord of the blade                                                     |
%                 (r) = radius station                                                         |
%                 (beta) = twist angle distribution                                            |
%                 (N) = number of blades                                                       |
%                 (First_lambda) = first value of tip speed                                    |
%                 (Last_lambda) = last value of tip speed                                      |
%                                                                                              |
% |Output      :  (Cp) = vector of power coefficient values for the                            |
%                 windmill                                                                     |
%                 (Cq) = vector of torque coefficient values for the                           |
%                 windmill                                                                     |
%                 (Ct) = vector of thrust coefficient values for the                           |
%                 windmill                                                                     |
%                 (Lambda) = vector of tip speed values corresponding to                       |
%                 the cp and cq coefficients                                                   |
% |Note        :                                                                               |
% ==============================================================================================


function [Cp, Cq, Ct, Lambda]=Curve_caratteristiche_aeromotori_asse_orizzontal(N,r,beta,c,First_lambda,Last_lambda)

R=r(end);
rs=r/R;

v_alpha=convang(linspace(1,60,100),'deg','rad');

cd=0.01;
Cl_alpha=6.28;

for i=1:length(r)
    for j=1:length(v_alpha)
        Cl(i,j)=Cl_alpha*v_alpha(j);

        phi(i,j)=v_alpha(j)+beta(i);
        Cx(i,j)=Cl(i,j)*cos(phi(i,j))+cd*sin(phi(i,j));
        Cy(i,j)=Cl(i,j)*sin(phi(i,j))-cd*cos(phi(i,j));
        sigmar(i,j)=N*c(i)/(2*pi*r(i));
        a(i,j)=(sigmar(i,j)*Cx(i,j)/(4*(sin(phi(i,j)))^2))/(1+(sigmar(i,j)*Cx(i,j)/(4*(sin(phi(i,j)))^2)));
        ap(i,j)=(sigmar(i,j)*Cy(i,j)/(4*sin(phi(i,j))*cos(phi(i,j))))/(1-(sigmar(i,j)*Cy(i,j)/(4*sin(phi(i,j))*cos(phi(i,j)))));
        lambda(i,j)=((1-a(i,j))/(1+ap(i,j)))*(1/tan(phi(i,j)))*(1/rs(i));

        dCtdrs(i,j)=2*sigmar(i,j)*Cx(i,j)*(r(i)/R)*(1-a(i,j))^2/((sin(phi(i,j)))^2);
        dCqdrs(i,j)=2*sigmar(i,j)*Cy(i,j)*(r(i)/R)^2*(1-a(i,j))^2/((sin(phi(i,j)))^2);
%        dCqdrs(i,j)=2*sigmar(i,j)*Cy(i,j)*(lambda(i,j))^2*(r(i)/R)^4*(1+ap(i,j))^2/((cos(phi(i,j)))^2);
        dCpdrs(i,j)=dCqdrs(i,j)*lambda(i,j);
        
    end

end



Lambda=linspace(First_lambda,Last_lambda,150);
STOP=ones(length(r),length(Lambda));
for k=1:length(Lambda)
    LAM=Lambda(k);
    for i=1:length(r)
        aa(i,k)=interp1(lambda(i,:),a(i,:),LAM,'pchip');
        aap(i,k)=interp1(lambda(i,:),ap(i,:),LAM,'pchip');
        
        if aa(i,k)>0.5
            STOP(i,k)=k;
        else
            STOP(i,k)=length(Lambda);
        end
            
        DCtdrs(i,k)=interp1(lambda(i,:),dCtdrs(i,:),LAM,'pchip');
        DCqdrs(i,k)=interp1(lambda(i,:),dCqdrs(i,:),LAM,'pchip');
        DCpdrs(i,k)=interp1(lambda(i,:),dCpdrs(i,:),LAM,'pchip');
 
    end
    
    Ct(k)=trapz(rs,DCtdrs(:,k));
    Cq(k)=trapz(rs,DCqdrs(:,k));
    Cp(k)=Cq(k)*LAM;
    
end

end


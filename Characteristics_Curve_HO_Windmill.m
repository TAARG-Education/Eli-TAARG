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
% |Name        : Characteristics_Curve_HO_Windmill.m                                           |
% |Author      : Nasti Giuseppe, Tomasso Armando Diego                                         |
% |              University of Naples Federico II.                                             |
% |Version     : 1.1                                                                           |
% |Date        : 20/01/2021                                                                    |
% |Modified    : 14/02/2021                                                                    |
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
%                 (Aeromatrix)= matrix of aerodynamic characteristics                          |                                                            |
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

function [Cp Cq Ct Lambda]=Characteristics_Curve_HO_Windmill(N,r,Aero_matrix,beta,c,firstlambda,lastlambda)

%% Input management

R=r(end); %Radius of the blade
rs=r/R; % adimensional radius station
v_alpha=convang(linspace(2,50,100),'deg','rad');
beta=convang(beta,'deg','rad'); %Blade Pitch 

%Extrapolation of the aerodynamic characteristics 
Cd_min    = Aero_matrix(:,1);
dCd_dCl2  = Aero_matrix(:,2);
Cl_Cd_min = Aero_matrix(:,3);
Re_ref    = Aero_matrix(:,4);
Re_inf    = Aero_matrix(:,5);
f         = Aero_matrix(:,6);
Cl_max    = Aero_matrix(:,7);
Cl_min    = Aero_matrix(:,8);
alpha_zero_lift=Aero_matrix(:,9); %      [deg]
cl_alpha=Aero_matrix(:,10); %            [1/deg]
dcldalfastall=Aero_matrix(:,11); %       [1/deg]

%% Gradient coefficents calculation

for i=1:length(r)
    
    %Creation of input vectors for the two subroutine lift_curve and ClCd_XRotor 
    v_input_cl(i,:)=[Cl_max(i),Cl_min(i),alpha_zero_lift(i),cl_alpha(i),dcldalfastall(i)];
    v_input_cd(i,:)=[Cd_min(i),dCd_dCl2(i),Cl_Cd_min(i),Re_ref(i),Re_inf(i),f(i),Cl_max(i),Cl_min(i)];
    
    for j=1:length(v_alpha)
        
        %Calculation of Cl and Cd for the blade station at angle v_alpha(i)
        Cl(i,j)=lift_curve(v_input_cl(i,:),convang(v_alpha(j),'rad','deg'));
        [~,~,Cd(i,j)]=ClCd_XRotor(v_input_cd(i,:), Cl(i,j));
        
        phi(i,j)=v_alpha(j)+beta(i); %inflow angle
        
        Cx(i,j)=Cl(i,j)*cos(phi(i,j))+Cd(i,j)*sin(phi(i,j));
        Cy(i,j)=Cl(i,j)*sin(phi(i,j))-Cd(i,j)*cos(phi(i,j));
        
        sigmar(i,j)=N*c(i)/(2*pi*r(i)); %blade solidity 
       
        a(i,j)=(sigmar(i,j)*Cx(i,j)/(4*(sin(phi(i,j)))^2))/(1+(sigmar(i,j)*Cx(i,j)/(4*(sin(phi(i,j)))^2))); 
        ap(i,j)=(sigmar(i,j)*Cy(i,j)/(4*sin(phi(i,j))*cos(phi(i,j))))/(1-(sigmar(i,j)*Cy(i,j)/(4*sin(phi(i,j))*cos(phi(i,j)))));
        
        lambda(i,j)=((1-a(i,j))/(1+ap(i,j)))*(1/tan(phi(i,j)))*(1/rs(i));

        dCtdrs(i,j)=2*sigmar(i,j)*Cx(i,j)*(r(i)/R)*(1-a(i,j))^2/((sin(phi(i,j)))^2);
        dCqdrs(i,j)=2*sigmar(i,j)*Cy(i,j)*(r(i)/R)^2*(1-a(i,j))^2/((sin(phi(i,j)))^2);
        dCpdrs(i,j)=dCqdrs(i,j)*lambda(i,j);
        
    end

end
%% Calculation of the characteristics curve
Lambda=linspace(firstlambda,lastlambda,150);

for k=1:length(Lambda)
    LAM=Lambda(k);
    for i=1:length(r)
        
        %Interpolation 
        aa(i,k)=interp1(lambda(i,:),a(i,:),LAM,'pchip');
        aap(i,k)=interp1(lambda(i,:),ap(i,:),LAM,'pchip');
            
        DCtdrs(i,k)=interp1(lambda(i,:),dCtdrs(i,:),LAM,'pchip');
        DCqdrs(i,k)=interp1(lambda(i,:),dCqdrs(i,:),LAM,'pchip');
        DCpdrs(i,k)=interp1(lambda(i,:),dCpdrs(i,:),LAM,'pchip');
 
    end
    
    if aa(:,k)<0.5 %Condition for Momentum theory validation
        
        %Coefficient calculation
        Ct(k)=trapz(rs,DCtdrs(:,k));
        Cq(k)=trapz(rs,DCqdrs(:,k));
        Cp(k)=Cq(k)*LAM;
    else 
        break %exit of the cycle for invalidation of the Momentum Theory (a>0.5)
    end
    
end

%Lambda selection in which momentum theory is valid
Lambda=Lambda(1:length(Cp));

%% Plot section

[y1,y2,y3]=Outputfunction('Characteristics_Curve_HO_Windmill',Lambda,Ct,Lambda,Cq,Lambda,Cp);
 
 
 end




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
% |Name        : turbine_Darrieus_tubo_flusso_multiplo.m                                       |
% |Author      : Luciano Roncioni - Camilla Scotto di Carlo                                    |
% |              University of Naples Federico II.                                             |
% |Version     : 1.5                                                                           |
% |Date        : 18/12/2020                                                                    |
% |Modified    : 27/12/2020                                                                    |
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

function [cp,cq,lambdav] = Multiple_Streamtube_Darrieus_Turbines_Theory(alphamax,c,R,cla,N,cd)
%% Initialization of vectors used in the code
phiv = linspace(0,360,100); %phi domain [deg]
phiv = deg2rad(phiv); % [rad]
lambdav = linspace(0.1,14,100); %tip speed
cpv = zeros(numel(lambdav),1); %cp vector
cqv = zeros(numel(lambdav),1); %cq vector
a = zeros(1,numel(phiv)); %induction
v_vinf = zeros(numel(phiv),numel(phiv)); %v/v_inf
alpha = zeros(numel(phiv),numel(phiv)); %alpha
ind = zeros(numel(phiv),numel(phiv)); %inductions vector
%% Values for entering the loop
cp = 0;
j = 0;
%% a, CP, CQ computing
while cp >= 0
    j = j +1;
    lambda = lambdav(j);
    for i = 1 : numel(phiv)
        phi = phiv(i);
        eq = @(a) ((1-a).*a) - (c/(4*R))* ( (lambda + (1-a).*sin(phi)).^2 +...
            (1-a).^2.* cos(phi).^2).*cla.*atan2(((1-a).*cos(phi)),(lambda + ...
            (1-a).*sin(phi))).*(cos(phi+(atan2(((1-a).*cos(phi)),(lambda +...
            (1-a).*sin(phi))))./cos(phi))); %anonymous function
        a(1,i) = (fzero(eq,0.2));%find zero of the previous function
    end
    
    ind(j,:) = a;
% Check on the reliability of the computed induction
% When ind is NaN, the cycle is left
    if sum(isnan(ind(j,:))) >= 1
        disp('Warning: from this lambda on, the results are not reliable,');
        disp('so the results are not reported');
        cp = -1;
        j = j-1;
    else
        
        v_vinf(j,:) = sqrt((lambda + (1-a).*sin(phiv)).^2 + (1-a).^2 .*cos(phiv).^2);
        alpha(j,:) = atan2(((1-a).*cos(phiv)),(lambda + (1-a).*sin(phiv)));
        cost_p = (N*c*lambda)./(4*pi*R);
        % computing Cp - numerical integration
        cp = cost_p.* trapz(phiv, (v_vinf(j,:).^2*cla.*alpha(j,:).*sin(alpha(j,:)).*...
            (1-(cd./(cla.*alpha(j,:))).*cot(alpha(j,:)))));
        cq = cp/lambda; % Cq value
        
        % Allocating values in corresponding vectors
        cpv(j,1) = cp;
        cqv(j,1) = cq;
        
        %Condition to exit the cycle when lambda
        % is the last value of the
        %vector
        if lambda == lambdav(end)
            cp = -1;
        end
    end
end
%Allocating vector of max values
max_v = zeros(numel(lambdav),1);
alphadeg = rad2deg(alpha); % [deg]

%Find max alpha for each row of the matrix alpha
for h = 1 : numel(phiv)
    max_v(h) = max(alphadeg(h,:));
end

counterlambdamin = 0; %intial value for lambda min index

%check on stall angle
for f = 1 : numel(phiv)
    if alphamax < max_v(f)
        counterlambdamin = counterlambdamin + 1;
    end
end

%cp, cq and lambdav are downsized according to conditions
% of stall and positive power value
cq = cqv(counterlambdamin:j,1);
cp = cpv(counterlambdamin:j,1);
lambdav = lambdav(1,counterlambdamin:j);
end
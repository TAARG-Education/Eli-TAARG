%% \dmst.m
%  \brief: VAWT performance analisys routine based on DMST theory
%  \author: Gabriele Lucci
%  \version: 1.00
%
% Eli-TAARG is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public
% License as published by the Free Software Foundation; either
% version 3 of the License, or (at your option) any later version.
%
% Eli-TAARG is developed by the TAARG Educational organization for
% educational purposes only.
% Theoretical and Applied Aerodynamic Research Group - University of Naples 
% Federico II.
%
% Eli-TAARG GitHub link: <https://github.com/TAARG-Education/Eli-TAARG>
%
% Eli-TAARG is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
% General Public License for more details.
% <http://www.gnu.org/licenses/>.
%
% =========================================================================
% |Name        : dmst.m                                                    
% |Author      : Gabriele Lucci                                            
% |              University of Naples Federico II.                         
% |Version     : 1.00                                                      
% |Date        : 02/05/2022                                                
% |Modified    : -                                                         
% |Description : H-Darrieus Vertical Axis Wind Turbines (VAWT) performance 
% |              assessment through Double-Multiple Streamtube theory,
% |              based on Ion Paraschivoiu (2002) approach.
% |Reference   : Tognaccini, R. (2020). "Lezioni di aerodinamica dell'ala 
% |              Rotante - A.A. 2019/2020". Università degli studi di 
% |              Napoli "Federico II".
% |              Paraschivoiu, I. (2002). "Wind Turbine Design With 
% |              Emphasis on Darrieus Concept". École Polytechnique de 
% |              Montréal.
% |              Sheldahl R. E. and Klimas P. C. (1981). "Aerodynamic 
% |              Characteristics of Seven Symmetrical Airfoil Sections 
% |              Through 180-Degree Angle of Attack for Use in Aerodynamic 
% |              Analysis of Vertical Axis Wind Turbines". Sandia National
% |              Laboratories.
% |Input       : (n_st)     = integer, number of streamtubes pairs;
% |              (B)        = integer, number of blades;
% |              (c)        = double, blade chord length [m];
% |              (R)        = double, rotor radius [m];
% |              (lambda)   = double, Tip Speed Ratio;
% |              (Vinf)     = double, wind speed [m/s];
% |              (aeroflag) = string, 'simple' | 'real'. Choose between 
% |                           simplified and realistic blade aerodynamic 
% |                           behaviour, assuming blade airfoil is a NACA 
% |                           0012.
% |Output      : (lambda_flag_us) = boolean, 1 if a > 0.5 somewhere 
% |                                 upwind, 0 otherwise;
% |              (lambda_flag_ds) = boolean, 1 if a > 0.5 somewhere 
% |                                 downwind, 0 otherwise;
% |              (lambda_eff_us)  = n_st-by-1 double array, upwind local
% |                                 TSR;
% |              (lambda_eff_ds)  = n_st-by-1 double array, downwind 
% |                                 local TSR;
% |              (Vratiosq_us)    = n_st-by-1 double array, upwind local 
% |                                 velocity ratio squared;
% |              (Vratiosq_ds)    = n_st-by-1 double array, downwind 
% |                                 local velocity ratio squared;
% |              (counter_us)     = n_st-by-1 integer array, upwind loop
% |                                 iteration counter;
% |              (counter_ds)     = n_st-by-1 integer array, downwind 
% |                                 loop iteration counter;
% |              (alpha_us)       = n_st-by-1 double array, upwind local
% |                                 angle of attack;
% |              (alpha_ds)       = n_st-by-1 double array, downwind
% |                                 local angle of attack;
% |              (Re_us)          = n_st-by-1 double array, upwind local
% |                                 Reynolds number;
% |              (Re_ds)          = n_st-by-1 double array, downwind 
% |                                 local Reynolds number;
% |              (Cn_us)          = n_st-by-1 double array, upwind local
% |                                 normal force coefficient;
% |              (Cn_ds)          = n_st-by-1 double array, downwind 
% |                                 local normal force coefficient;
% |              (Ct_us)          = n_st-by-1 double array, upwind local
% |                                 tangential force coefficient;
% |              (Ct_ds)          = n_st-by-1 double array, downwind
% |                                 local tangential force coefficient;
% |              (a_us)           = n_st-by-1 double array, upwind local
% |                                 axial induction factor;
% |              (a_ds)           = n_st-by-1 double array, downwind
% |                                 local axial induction factor;
% |              (instCq_us)      = n_st-by-1 double array, instantaneous 
% |                                 torque coefficient for the upwind 
% |                                 cycle;
% |              (instCq_ds)      = n_st-by-1 double array, instantaneous 
% |                                 torque coefficient for the downwind 
% |                                 cycle;
% |              (CP_us)          = double, power coefficient generated by 
% |                                 B blades in the upwind passage, 
% |                                 averaged on the whole rotor revolution;
% |              (CP_ds)          = double, power coefficient generated by 
% |                                 B blades in the downwind passage, 
% |                                 averaged on the whole rotor revolution;
% |              (CP)             = double, average net power coefficient.
% |Note        : This function capabilities are limited to the aerodynamic
% |              data available (if the realistic aerodynamic model is
% |              chosen), and to the extent that the hypothesis underneath 
% |              Double-Multiple Streamtube theory are reasonably
% |              acceptable.
% =========================================================================
%
%%
function [ ...
    lambda_flag_us,lambda_flag_ds, ...
    lambda_eff_us,lambda_eff_ds, ...
    Vratiosq_us,Vratiosq_ds, ...
    counter_us,counter_ds, ...
    alpha_us,alpha_ds, ...
    Re_us,Re_ds, ...
    Cn_us,Cn_ds, ...
    Ct_us,Ct_ds, ...
    a_us,a_ds, ...
    instCQ_us,instCQ_ds, ...
    CP_us,CP_ds,CP ...
    ] = dmst(n_st,B,c,R,lambda,Vinf,aeroflag)

global RADrunflag RE

% Check wether realistic model has been chosen, but aerodynamic data was
% not previoulsy loaded.
if strcmpi(aeroflag,'real') && isempty(RADrunflag)
    
    filename = 'sandia0012data.xlsx';
    filepath = [cd,'/ExperimentalData/',filename];
    ReadAeroData(filepath);
    
end    

Delta_theta = pi/n_st;

theta_us_seq = linspace(pi/2,-pi/2,n_st);
theta_ds_seq = linspace(pi/2,3*pi/2,n_st);

sigma = B*c/R;
nu = 1.5e-5;

%% Vars init
lambda_eff_us = zeros(n_st,1);
lambda_eff_ds = zeros(n_st,1);
Vratiosq_us   = zeros(n_st,1);
Vratiosq_ds   = zeros(n_st,1);
counter_us    = zeros(n_st,1);
counter_ds    = zeros(n_st,1);
instCQ_us     = zeros(n_st,1);
instCQ_ds     = zeros(n_st,1);
alpha_us      = zeros(n_st,1);
alpha_ds      = zeros(n_st,1);
Re_us         = zeros(n_st,1);
Re_ds         = zeros(n_st,1);
Cn_us         = zeros(n_st,1);
Cn_ds         = zeros(n_st,1);
Ct_us         = zeros(n_st,1);
Ct_ds         = zeros(n_st,1);
u_us          = zeros(n_st,1);
u_ds          = zeros(n_st,1);

%% Calc loops
disp('Entering upwind loop...');
for ind_theta = 1:n_st
    
    theta = theta_us_seq(ind_theta);
    
    u_us_old = 1;
    
    exitflag = -1;
    
    while exitflag == -1
        
        counter_us(ind_theta) = counter_us(ind_theta) + 1;
        
        lambda_eff_us(ind_theta) = lambda/u_us_old;
        
        Vratiosq_us(ind_theta) = ...
            (lambda_eff_us(ind_theta) - ...
            sin(theta))^2 + cos(theta)^2;
        
        alpha_us(ind_theta) = ...
            asin(cos(theta)/ ...
            sqrt(Vratiosq_us(ind_theta)));
        
        Re_us(ind_theta) = ...
            Vinf*sqrt(Vratiosq_us(ind_theta))*c/nu;
        
        % Check if local Reynolds number outranges tabulated values
        if strcmpi(aeroflag,'real')
            
            if Re_us(ind_theta) < RE(1)
                
                error(['Local Reynolds number is lower than minimum ', ...
                    'value (',num2str(RE(1)),') in the available data.' ...
                    ' Cannot lookup aerodynamics table and continue.', ...
                    newline, 'theta = ',num2str(theta), ...
                    ', lambda = ',num2str(lambda), ...
                    ', Re = ',num2str(Re_us(ind_theta)),'.']);
                
            elseif Re_us(ind_theta) > RE(end)
                
                error(['Local Reynolds number is grater than maximum', ...
                    'value (',num2str(RE(end)), ...
                    ') in the available data. Cannot lookup  ',...
                    'aerodynamics table and continue.', ...
                    newline, 'theta = ',num2str(theta), ...
                    ', lambda = ',num2str(lambda), ...
                    ', Re = ',num2str(Re_us(ind_theta)),'.']);
                
            end
            
        end
        
        Cn_us(ind_theta) = ...
            Cl_dmst(Re_us(ind_theta), ...
            alpha_us(ind_theta), ...
            aeroflag)* ...
            cos(alpha_us(ind_theta)) + ...
            Cd_dmst(Re_us(ind_theta), ...
            alpha_us(ind_theta), ...
            aeroflag)* ...
            sin(alpha_us(ind_theta));
        
        Ct_us(ind_theta) = ...
            Cl_dmst(Re_us(ind_theta), ...
            alpha_us(ind_theta), ...
            aeroflag)* ...
            sin(alpha_us(ind_theta)) - ...
            Cd_dmst(Re_us(ind_theta), ...
            alpha_us(ind_theta), ...
            aeroflag)* ...
            cos(alpha_us(ind_theta));
        
        
        if ind_theta == 1 || ind_theta == n_st
            
            % Exit loop where F_us would be singular...
            exitflag = 1;
            u_us(ind_theta) = 1;
            
        else
            
            % ...or find new induction value
            intfun = @(theta) ...
                Vratiosq_us(ind_theta).* ...
                (Cn_us(ind_theta).*cos(theta)./ ...
                abs(cos(theta)) - ...
                Ct_us(ind_theta).*sin(theta)./ ...
                abs(cos(theta)));
            
            F_us = sigma/ ...
                (8*Delta_theta)* ...
                integral(intfun, ...
                theta-Delta_theta/2,theta+Delta_theta/2);
            
            u_us_new = pi/(F_us + pi);
            
            % Convergence check
            if abs(u_us_old - u_us_new) < 1e-2
                
                exitflag = 1;
                
                % Update variables
                u_us(ind_theta) = u_us_new;
                
                lambda_eff_us(ind_theta) = lambda/u_us_new;
                
                Vratiosq_us(ind_theta) = ...
                    (lambda_eff_us(ind_theta) - ...
                    sin(theta))^2 + cos(theta)^2;
                
                alpha_us(ind_theta) = ...
                    asin(cos(theta)/ ...
                    sqrt(Vratiosq_us(ind_theta)));
                
                Re_us(ind_theta) =  ...
                    Vinf*sqrt(Vratiosq_us(ind_theta))*c/nu;
                
                Cn_us(ind_theta) = ...
                    Cl_dmst(Re_us(ind_theta), ...
                    alpha_us(ind_theta), ...
                    aeroflag)* ...
                    cos(alpha_us(ind_theta)) + ...
                    Cd_dmst(Re_us(ind_theta), ...
                    alpha_us(ind_theta), ...
                    aeroflag)* ...
                    sin(alpha_us(ind_theta));
                
                Ct_us(ind_theta) =  ...
                    Cl_dmst(Re_us(ind_theta), ...
                    alpha_us(ind_theta), ...
                    aeroflag)* ...
                    sin(alpha_us(ind_theta)) - ...
                    Cd_dmst(Re_us(ind_theta), ...
                    alpha_us(ind_theta), ...
                    aeroflag)* ...
                    cos(alpha_us(ind_theta));
                
            end
            
            u_us_old = u_us_new;
            
        end
        
    end
    
    instCQ_us(ind_theta) = ...
        sigma/4*u_us(ind_theta)^2*Vratiosq_us(ind_theta)*Ct_us(ind_theta);
    
end

CQ_us = sigma/(8*pi)* ...
    trapz(flip(theta_us_seq), ...
    Ct_us(:).*u_us(:).^2.* ...
    Vratiosq_us(:));

CP_us = CQ_us*lambda;

a_us = 1 - u_us;

disp(['...upwind problem solved in ', ...
    num2str(sum(counter_us)),' total iterations.']);

if any(a_us > 0.5)
    
    lambda_flag_us = 1;
    warning(['Rankine-Froude theory limit (a = 0.5) exceeded ', ...
        '<strong>upwind</strong>.',newline, ...
        'Entering downwind loop...']);
    
else
    
    lambda_flag_us = 0;
    disp('Entering downwind loop...');
    
end

for ind_theta = 1:n_st
    
    theta = theta_ds_seq(ind_theta);
    
    u_ds_old = u_us(ind_theta);
    
    exitflag = -1;
    
    while exitflag == -1
        
        counter_ds(ind_theta) = counter_us(ind_theta) + 1;
        
        u_us_local = u_us(ind_theta);
        
        lambda_eff_ds(ind_theta) = lambda/(2*u_us_local - 1)*u_ds_old;
        
        Vratiosq_ds(ind_theta) = ...
            (lambda_eff_ds(ind_theta) - ...
            sin(theta))^2 + cos(theta)^2;
        
        alpha_ds(ind_theta) = ...
            asin(cos(theta)/ ...
            sqrt(Vratiosq_ds(ind_theta)));
        
        Re_ds(ind_theta) = ...
            Vinf*sqrt(Vratiosq_ds(ind_theta))/nu;
        
        % Check if local Reynolds number outranges tabulated values
        if strcmpi(aeroflag,'real')
            
            if Re_ds(ind_theta) < RE(1)
                
                error(['Local Reynolds number is lower than minimum ', ...
                    'value (',num2str(RE(1)),') in the available data.' ...
                    ' Cannot lookup aerodynamics table and continue.', ...
                    newline, 'theta = ',num2str(theta), ...
                    ', lambda = ',num2str(lambda), ...
                    ', Re = ',num2str(Re_ds(ind_theta)),'.']);
                
            elseif Re_ds(ind_theta) > RE(end)
                
                error(['Local Reynolds number is grater than maximum', ...
                    'value (',num2str(RE(end)), ...
                    ') in the available data. Cannot lookup  ',...
                    'aerodynamics table and continue.', ...
                    newline, 'theta = ',num2str(theta), ...
                    ', lambda = ',num2str(lambda), ...
                    ', Re = ',num2str(Re_ds(ind_theta)),'.']);
                
            end
            
        end
        
        Cn_ds(ind_theta) = ...
            Cl_dmst(Re_ds(ind_theta), ...
            alpha_ds(ind_theta), ...
            aeroflag)* ...
            cos(alpha_ds(ind_theta)) + ...
            Cd_dmst(Re_ds(ind_theta), ...
            alpha_ds(ind_theta), ...
            aeroflag)* ...
            sin(alpha_ds(ind_theta));
        
        Ct_ds(ind_theta) = ...
            Cl_dmst(Re_ds(ind_theta), ...
            alpha_ds(ind_theta), ...
            aeroflag)* ...
            sin(alpha_ds(ind_theta)) - ...
            Cd_dmst(Re_ds(ind_theta), ...
            alpha_ds(ind_theta), ...
            aeroflag)* ...
            cos(alpha_ds(ind_theta));
        
        if ind_theta == 1 || ind_theta == n_st
            
            % Exit loop where F_ds would be singular
            exitflag = 1;
            u_ds(ind_theta) = 1;
            
        else
            
            % ...or find new induction value
            intfun = @(theta) ...
                Vratiosq_ds(ind_theta).* ...
                (Cn_ds(ind_theta).*cos(theta)./ ...
                abs(cos(theta)) - ...
                Ct_ds(ind_theta).*sin(theta)./ ...
                abs(cos(theta)));
            
            F_ds = sigma/(8*Delta_theta)* ...
                integral(intfun, ...
                theta-Delta_theta/2,theta+Delta_theta/2);
            
            u_ds_new = pi/(F_ds + pi);
            
            % Convergence check
            if abs(u_ds_old - u_ds_new) < 1e-2
                
                exitflag = 1;
                
                % Update variables
                u_ds(ind_theta) = u_ds_new;
                
                lambda_eff_ds(ind_theta) = ...
                    lambda/(2*u_us_local - 1)*u_ds_new;
                
                Vratiosq_ds(ind_theta) = ...
                    (lambda_eff_ds(ind_theta) - ...
                    sin(theta))^2 + cos(theta)^2;
                
                alpha_ds(ind_theta) = ...
                    asin(cos(theta)/ ...
                    sqrt(Vratiosq_ds(ind_theta)));
                
                Re_ds(ind_theta) = ...
                    Vinf*sqrt(Vratiosq_ds(ind_theta))*c/nu;
                
                Cn_ds(ind_theta) = ...
                    Cl_dmst(Re_ds(ind_theta), ...
                    alpha_ds(ind_theta), ...
                    aeroflag)* ...
                    cos(alpha_ds(ind_theta)) + ...
                    Cd_dmst(Re_ds(ind_theta), ...
                    alpha_ds(ind_theta), ...
                    aeroflag)* ...
                    sin(alpha_ds(ind_theta));
                
                Ct_ds(ind_theta) = ...
                    Cl_dmst(Re_ds(ind_theta), ...
                    alpha_ds(ind_theta), ...
                    aeroflag)* ...
                    sin(alpha_ds(ind_theta)) - ...
                    Cd_dmst(Re_ds(ind_theta), ...
                    alpha_ds(ind_theta), ...
                    aeroflag)* ...
                    cos(alpha_ds(ind_theta));
                
            end
            
            u_ds_old = u_ds_new;
            
        end
        
    end
    
    instCQ_ds(ind_theta) = ...
        sigma/4*...
        (u_ds(ind_theta)*(2*u_us(ind_theta) - 1))^2* ...
        Vratiosq_ds(ind_theta)*Ct_ds(ind_theta);    
    
end

CQ_ds = sigma/(8*pi)*trapz(theta_ds_seq, ...
    Ct_ds(:).* ...
    ((2*u_us(ind_theta) - 1)*u_ds(:)).^2.* ...
    Vratiosq_ds(:));

CP_ds = CQ_ds*lambda;

CP = CP_us + CP_ds;

a_ds = 1 - u_ds;

disp(['...downwind problem solved in ', ...
    num2str(sum(counter_ds)),' total iterartions.']);

if any(a_ds > 0.5)
    
    lambda_flag_ds = 1;
    warning(['Rankine-Froude theory limit (a = 0.5) exceeded ', ...
        '<strong>downwind</strong>.']);
    
else
    
    lambda_flag_ds = 0;
    
end

end
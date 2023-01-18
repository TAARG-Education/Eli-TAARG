%% Helicopter SIKORSKY UH-60A %%

 clc; close all; clear all;

 % DATA


  % 1) aircraft
    Wto        = 9979.03*9.81;     % [N]
    Pd_sl        = 2110000;        % [W] 
    h          = 0;                % [m]

    T0 = 288; % K
    mu0 = 1.79e-5; % Kg/(ms)
    [Temp, a0, P0, rho] = atmosisa(h);
    mu_viscoso = mu0*(Temp/T0)^1.5*(T0+110)/(Temp+110);  % Sutherlan's law
 
  % 2) main rotor                      % tail rotor
    N        = 4;                       N_rc     = 4;
    R_rp     = 8.18;                    R_rc     = 1.7;          % [m]
    A_rp     = pi*R_rp^2;               A_rc     = pi*R_rc^2;    % [m^2]
    c        = 0.53;                    c_rc     = 0.2;          % [m]
    sigma_rp = 0.082;                   sigma_rc = 0.188;
    Vtip     = 220.98;                  Vtip_rc  = 208.79;       % [m/s]
    Omega_rp = Vtip/R_rp;               Omega_rc = Vtip_rc/R_rc;
    n        = Omega_rp/(2*pi);         n_rc     = Omega_rc/(2*pi);
  
                           br = 10.73;       % [m]
    
    theta1_rp = -18;                      theta1_rc = -17;         % [deg]
    Num       = 10;
    r_rp      = linspace(0.05,R_rp,Num);  r_rc      = linspace(0.05,R_rc,Num);
    r_adim_rp = r_rp/R_rp;                r_adim_rc = r_rc./R_rc;

  Cd_rc=0.01;
  Cd_rp=0.01;
 
 k = 1.15;                  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Forward flight - request power at differente altitudes and weights %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


numero  = 120;
Vinf    = linspace(0,120,numero);
P_aus   = 25000;  % [W]
P_ausil = P_aus*ones(1,length(Vinf));
eta_rc   = 1.03;   

%available power  
v_h = linspace(0, 1800, length(Vinf));
[~, ~, ~, rho_sl] = atmosisa(0);
for j = 1:length(v_h)
  [a,b,c,density] = atmosisa(v_h(j));
  v_Pdisponibile(j)= Pd_sl*density/rho_sl;
end
figure 
plot(v_h,v_Pdisponibile, 'k'); 
xlabel('h [m]'); ylabel('P_a_v_a_i [W]'); title('Avaiable Power at different h');
grid on;


%%% total power at different altitudes and weights
v_WTO = [0.75*Wto,0.90*Wto, Wto];       
h = [0,1000,1800];



v_input=[R_rp,A_rp,sigma_rp, Omega_rp, Cd_rp,...
       R_rc, A_rc, sigma_rc, Omega_rc, Cd_rc,...
         br, k, P_aus, eta_rc, Pd_sl ];

for i = 1:length(h)
    for j = 1:length(v_WTO)
    
    [Pindotta_B,  Pparassita_rp,Pfus,Ptot_rp,...
     Pindotta_rc,  Pparassita_rc, Ptot_rc,Ptot_richiesta(j,:),Pdisp,...
     V_NE,V_BE, V_BR] = PowerInForwardFlight_main_tail_rotor(h(i), Vinf, v_WTO(j),v_input);
    
    Pdisponibile(j,:) = Pdisp*ones(1,length(Vinf));
 
    if (h(i) == 0   &&   v_WTO(j) == Wto)
   
    figure
    plot(Vinf,Pindotta_B,'k-.',Vinf,Pparassita_rp,'k--',Vinf,Pfus,'k:',Vinf,Ptot_rp,'k-',Vinf,Pdisponibile(j,:),'.-k');
    hold on;

    xlabel('V_∞ [m/s]'); ylabel('P [W]'); title('Required Power of the main rotor in forward flight');
    subtitle(['altitude h = ',num2str(h(i)), 'm']);
    legend('P_i_n_d','P_p_a_r','P_f_u_s_o_l_a_g_e','P_t_o_t','P_a_v_a_i');
    grid on
    


    figure
    plot(Vinf,Pparassita_rc,'k-.',Vinf,Pindotta_rc,'k--',Vinf,Ptot_rc,'k-');
    title('Powers of the tail rotor'); xlabel('V_∞ [m/s]'); ylabel('P [W]');
    subtitle(['altitude h = ',num2str(h(i)), 'm']);
    legend('P_p_a_r','P_i_n_d','P_t_o_t');
     grid on;

     figure
     plot(Vinf,Ptot_rp,'k-.',Vinf,Ptot_rc,'k--',Vinf,P_ausil,'k:',Vinf,Ptot_richiesta(j,:),'k-',Vinf,Pdisponibile(j,:),'.-k');
     hold on;
     xlabel('V_∞ [m/s]'); ylabel('P [W]'); title('Total Required Power in forward flight');
     subtitle(['altitude h = ',num2str(h(i)), 'm']);
     legend('P_m_a_i_n _r_o_t_o_r.','P_t_a_i_l _r_o_t_o_r','P_a_u_x','P_t_o_t','P_a_v_a_i');
     grid on
 
    end 
 end   

   figure
    plot(Vinf,Ptot_richiesta(j,:),'.-k');
    xlabel('V_∞ [m/s]'); ylabel('P [W]');
    subtitle(['altitude h = ',num2str(h(i)), 'm']);
    hold on;
    plot(Vinf,Pdisponibile(j,:),'k');
    assey = linspace(0,Ptot_richiesta(end),numero);
    plot(V_NE*ones(1,numero),assey,'k-.');
    plot(V_BE*ones(1,numero),assey,'k--');
    plot(V_BR*ones(1,numero),assey,'k');
    legend('P_t_o_t','P_a_v_a_i','V_N_E','V_B_E','V_B_R');

 % Total required power  in forward flight at different WTO
    figure
      plot(Vinf, Ptot_richiesta(1,:), 'k- ',...
      Vinf, Ptot_richiesta(2,:), 'k--',...
      Vinf, Ptot_richiesta(3,:), 'k-.',...
      Vinf, Pdisponibile(j,:)  , '.-k ');

    title('Total required power in forward flight at different W_T_O')
    subtitle(['altitude h = ',num2str(h(i)), 'm']);
    xlabel('V_∞ [m/s]'); ylabel('P [W]');
    legend('75% WTO',...
           '90% WTO',...
            'WTO',...
            'P_a_v_a_i_l_a_b_l_e');

    grid on;

end



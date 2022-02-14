close all; clear all; clc;

%% Setup
n_st = 25;
theta_us_seq = linspace(pi/2,-pi/2,n_st);
theta_ds_seq = linspace(pi/2,3*pi/2,n_st);
theta_seq    = [fliplr(theta_us_seq),theta_ds_seq(2:end)];

lambda_seq = 1.5:0.05:5.8;
% lambda_seq = 1:0.05:6;
n_lambda   = length(lambda_seq);

% VAWT rotor geometry
c = 0.2;
R = 2;
B = 3;
sigma = B*c/R;

% Aerodynamics
aeroflag = 'real';
Vinf = 5;

%% Vars init
lambda_flag_us = zeros(n_lambda,1);
lambda_flag_ds = zeros(n_lambda,1);
lambda_eff_us  = zeros(n_st,n_lambda);
lambda_eff_ds  = zeros(n_st,n_lambda);
Vratiosq_us    = zeros(n_st,n_lambda);
Vratiosq_ds    = zeros(n_st,n_lambda);
counter_us     = zeros(n_st,n_lambda);
counter_ds     = zeros(n_st,n_lambda);
instCQ_us      = zeros(n_st,n_lambda);
instCQ_ds      = zeros(n_st,n_lambda);
alpha_us       = zeros(n_st,n_lambda);
alpha_ds       = zeros(n_st,n_lambda);
Re_us          = zeros(n_st,n_lambda);
Re_ds          = zeros(n_st,n_lambda);
Cn_us          = zeros(n_st,n_lambda);
Cn_ds          = zeros(n_st,n_lambda);
Ct_us          = zeros(n_st,n_lambda);
Ct_ds          = zeros(n_st,n_lambda);
CQ_us          = zeros(n_lambda,1);
CQ_ds          = zeros(n_lambda,1);
CP_us          = zeros(n_lambda,1);
CP_ds          = zeros(n_lambda,1);
a_us           = zeros(n_st,n_lambda);
a_ds           = zeros(n_st,n_lambda);
CP             = zeros(n_lambda,1);

%% Performance calculation at various TSR
disp(['<strong>VAWT performance analysis ', ...
    'through Double-Multiple Streamtube method</strong>',newline, ...
    'sigma = ',num2str(sigma),', Vinf = ',num2str(Vinf),' m/s']);

if strcmpi(aeroflag,'simple')
    
    aeromessage = ...
        ['Simplified aerodynamic model.',newline];
    
elseif strcmpi(aeroflag,'real')
    
    aeromessage = ...
        ['Aerodynamics from NACA 0012 Sheldahl & Klimas data.',newline];
    
else
    
    error("Spellcheck 'aeroflag'")
    
end

disp(aeromessage);

for ind_lambda = 1:n_lambda
    
    lambda = lambda_seq(ind_lambda);
    
    disp(['Solving for lambda = ',num2str(lambda)'.']);
    
    [ ...
        lambda_flag_us(ind_lambda), ...
        lambda_flag_ds(ind_lambda), ...
        lambda_eff_us(:,ind_lambda), ...
        lambda_eff_ds(:,ind_lambda), ...
        Vratiosq_us(:,ind_lambda), ...
        Vratiosq_ds(:,ind_lambda), ...
        counter_us(:,ind_lambda), ...
        counter_ds(:,ind_lambda), ...
        alpha_us(:,ind_lambda), ...
        alpha_ds(:,ind_lambda), ...
        Re_us(:,ind_lambda), ...
        Re_ds(:,ind_lambda), ...
        Cn_us(:,ind_lambda), ...
        Cn_ds(:,ind_lambda), ...
        Ct_us(:,ind_lambda), ...
        Ct_ds(:,ind_lambda), ...
        a_us(:,ind_lambda), ...
        a_ds(:,ind_lambda), ...
        instCQ_us(:,ind_lambda),instCQ_ds(:,ind_lambda), ...
        CP_us(ind_lambda),CP_ds(ind_lambda),CP(ind_lambda) ...
        ] = dmst(n_st,B,c,R,lambda,Vinf,aeroflag);
    
    disp([newline,'--------------------------------------------------', ...
        newline]);
    
end

%% Diagrams
% Power coefficient
f1 = figure(1);
hold on;

if any(lambda_flag_us) || any(lambda_flag_ds)
    
    ind_lambda_flag_us = find(lambda_flag_us,1);
    ind_lambda_flag_ds = find(lambda_flag_ds,1);
    ind_lambda_flag = min([ind_lambda_flag_us,ind_lambda_flag_ds]);
    plot(lambda_seq(1:ind_lambda_flag),CP(1:ind_lambda_flag),'k');
    invalid_vals = ...
        plot(lambda_seq(ind_lambda_flag:end), ...
        CP(ind_lambda_flag:end),'--k');
    invalid_vals.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
else
    
    plot(lambda_seq,CP,'k');
    
end

title('Double-Multiple Streamtube theory - C_P vs. \lambda', ...
    ['Dashed line: local MT breakdown (a > 0.5).', ...
    ' Rotor \sigma = ',num2str(sigma)]);
xlabel('\lambda'); ylabel('C_P'); grid minor;
xlim([lambda_seq(1),lambda_seq(end)]);

% Power coefficient breakdown
f2 = figure(2);
hold on;

plot(lambda_seq,CP_us,'-.k',lambda_seq,CP_ds,'--k',lambda_seq,CP,'k');

legend('C_{P,us}','C_{P,ds}','C_{P,net}','Location','northwest');
title('Power coefficient breakdown')
xlabel('\lambda'); ylabel('C_P'); grid minor;
xlim([lambda_seq(1),lambda_seq(end)]);

% Axial inductions
f3 = figure(3);
a = [flipud(a_us);a_ds(2:end,:)];

plot(theta_seq,a(:,1),'-.k', ...
    theta_seq,a(:,floor(n_lambda/2)),'--k', ...
    theta_seq,a(:,end),'k'); hold on;
plot(theta_seq,0.5*ones(size(theta_seq)),'k');

legend(['\lambda = ',num2str(lambda_seq(1))], ...
    ['\lambda = ',num2str(lambda_seq(floor(n_lambda/2)))], ...
    ['\lambda = ',num2str(lambda_seq(end))],'', ...
    'Location','northeastoutside');
title(['Axial induction, \sigma = ',num2str(sigma),'.']);
xlabel('\theta'); ylabel('a'); axis tight; grid minor; 

% Angle of attack
f4 = figure(4);
alpha = [flipud(alpha_us);alpha_ds(2:end,:)];

plot(theta_seq,alpha(:,1),'-.k', ...
    theta_seq,alpha(:,floor(n_lambda/2)),'--k', ...
    theta_seq,alpha(:,end),'k'); hold on;

legend(['\lambda = ',num2str(lambda_seq(1))], ...
    ['\lambda = ',num2str(lambda_seq(floor(n_lambda/2)))], ...
    ['\lambda = ',num2str(lambda_seq(end))], ...
    'Location','northeastoutside');
title(['Angle of attack, \sigma = ',num2str(sigma),'.']);
xlabel('\theta'); ylabel('\alpha'); axis tight; grid minor; 

% Instantaneous torque coefficient
f5 = figure(5);
instCQ = [flipud(instCQ_us);instCQ_ds(2:end,:)];

plot(theta_seq,instCQ(:,1),'-.k', ...
    theta_seq,instCQ(:,61),'k'); hold on;

legend(['\lambda = ',num2str(lambda_seq(1))], ...
    ['\lambda = ',num2str(lambda_seq(61))], ...
    'Location','northeastoutside');
title(['Instantaneous torque coefficient, \sigma = ',num2str(sigma),'.']);
xlabel('\theta'); ylabel('C_Q'); axis tight; grid minor; 
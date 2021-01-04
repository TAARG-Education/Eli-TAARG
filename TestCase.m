%This script is used to test the function. The value are compared with the
%one obtained from: R. Tognaccini. "Lezioni di aerodinamica dell’ala rotante" 2019 pp. 84-85,fig. 7.1-4. 

clear all;close all;clc;

%% Import data from text file
%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["V", "w"];
opts.VariableTypes = ["double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
tbl = readtable("C:\Users\dtagl\Desktop\TestCase\InductionThrustCalfa0.txt", opts);

%% Convert to output type
Vt_t1 = tbl.V;
wt_t = tbl.w;

%% Clear temporary variables
clear opts tbl

%% Import data from text file


%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["Vt", "Pt"];
opts.VariableTypes = ["double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
tbl = readtable("C:\Users\dtagl\Desktop\TestCase\PowerThrustCalfa0.txt", opts);

%% Convert to output type
Vt_t2 = tbl.Vt;
Pt_t = tbl.Pt;

%% Clear temporary variables
clear opts tbl

%% Import data from text file

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["Vp", "wp"];
opts.VariableTypes = ["double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
tbl = readtable("C:\Users\dtagl\Desktop\TestCase\InductionPowerCalfa0.txt", opts);

%% Convert to output type
Vp_t1 = tbl.Vp;
wp_t = tbl.wp;

%% Clear temporary variables
clear opts tbl

%% Import data from text file

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["Vp", "Tp"];
opts.VariableTypes = ["double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
tbl = readtable("C:\Users\dtagl\Desktop\TestCase\ThrustPowerCalfa0.txt", opts);

%% Convert to output type
Vp_t2 = tbl.Vp;
Tp_t = tbl.Tp;

%% Clear temporary variables
clear opts tbl

[Vt,wt,Vp,wp,Pt,Tp] = RotorFF(0);


subplot(2,2,1);
plot(Vt,wt,'k',Vt_t1,wt_t,'k--');
axis square;
xlabel('$\tilde{V_{\infty}}$','Interpreter','latex');ylabel('$\tilde{w}$','Interpreter','latex');
title('TestCase curve $\tilde{w} = \tilde{w}(\tilde{V_{\infty}})$,constant Thrust', 'Interpreter','latex') ;
legend('Function','Test Reference');


subplot(2,2,3);
plot(Vt,Pt,'k',Vt_t2,Pt_t,'k--');
axis square;
xlabel('$\tilde{V_{\infty}}$','Interpreter','latex');ylabel('$\tilde{P}$','Interpreter','latex');
title('TestCase curve $\tilde{P} = \tilde{P}(\tilde{V_{\infty}})$,constant Thrust', 'Interpreter','latex') ;
legend('Function','Test Reference');


subplot(2,2,2);
plot(Vp,wp,'k',Vp_t1,wp_t,'k--');
axis square;axis([0 5 0 1]);
xlabel('$\tilde{V_{\infty}}$','Interpreter','latex');ylabel('$\tilde{w}$','Interpreter','latex');
title('TestCase curve $\tilde{w} = \tilde{w}(\tilde{V_{\infty}})$,constant Power', 'Interpreter','latex') ;
legend('Function','Test Reference');


subplot(2,2,4);
plot(Vp,Tp,'k',Vp_t2,Tp_t,'k--');
axis square;axis([0 5 0 3]);
xlabel('$\tilde{V_{\infty}}$','Interpreter','latex');ylabel('$\tilde{T}$','Interpreter','latex');
title('TestCase curve $\tilde{T} = \tilde{T}(\tilde{V_{\infty}})$,constant Power', 'Interpreter','latex') ;
legend('Function','Test Reference');
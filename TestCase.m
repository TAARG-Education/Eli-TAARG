%% \RotorFF
%  \The function plots the characteristic curves for rotor in forward
%  flight for both constant thrust and power and gives in output also the
%  relative x and y values. It requires in input the angle of attack in 
%  degrees.
%  The plot avaiable are:
%- For constant Thrust:
%  -w versus V
%  -P versus V
%- For constant Power:
%  -w versus V
%  -T versus V
% where w = induction, V = asymptotic velocity, T = Thrust, P = Power
% All the values are non-dimensional with respect to their value in
% hovering (for V is used induction in hovering).
%  \Angelo Tagliafierro
%  \0.1
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
% |Name        : RotorFF                                                                       |
% |Author      : Angelo Tagliafierro                                                           |
% |              University of Naples Federico II.                                             |
% |Version     : 0.1                                                                           |
% |Date        : 15/12/2020                                                                    |
% |Modified    : 15/12/2020                                                                    |
% |Description : The function generates the characteristic curves for rotor in forward         |
%                flight and the relative x and y values                                        |
% |Reference   : R. Tognaccini. "Lezioni di aerodinamica dell’ala rotante" 2019 pp. 84-85.     |
% |Input       : Angle of attack                                                               |
% |Output      : Plots and axis values                                                         |
% |Note        :                                                                               |
% ==============================================================================================

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
tbl = readtable("InductionThrustCalfa0.txt", opts);

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
tbl = readtable("PowerThrustCalfa0.txt", opts);

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
tbl = readtable("InductionPowerCalfa0.txt", opts);

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
tbl = readtable("ThrustPowerCalfa0.txt", opts);

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

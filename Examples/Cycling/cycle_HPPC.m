% -------------------------------------------------------------------------
% cycle_HPPC.m
% -------------------------------------------------------------------------
%   cycle_HPPC is an example code showing how to simulate HPPC profile. 
% 
%   HPPC consists of simulating current pulses and rest periods during CC discharging.
%   Here we simulate HPPC by concatenating multiple CC simulations, called CC segments.
%   For more information on HPPC tests, refer to:
%   [1] S. Ha, G. Pozzato, and S. Onori, "Electrochemical characterization tools for lithium-ion batteries," J Solid State Electrochem, Nov. 2023, doi: 10.1007/s10008-023-05717-1.
%   [2] G. Pozzato, A. Allam, and S. Onori, "Lithium-ion battery aging dataset based on electric vehicle real-driving profiles," Data in Brief, vol. 41, p. 107995, Apr. 2022, doi: 10.1016/j.dib.2022.107995.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COBRAPRO: Co-simulation Battery Modeling for Accelerated Parameter Optimization

% Copyright (c) 2024 CO-simulation BatteRy modeling for Accelerated PaRameter Optimization (COBRAPRO)
% COBRAPRO is freely distributed under the MIT License 
% v1.0.0.: Released March, 2024 

% Main contributors: 
% Sara Ha (sungyeon.sara.ha@stanford.edu)
% Simona Onori (sonori@stanford.edu)
% Stanford Energy Control Group, Stanford University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% User Input  
clc;clear;close all

%--------------------------------------------------------------------------
% Option 1: Enter identified parameters
%--------------------------------------------------------------------------
load('identified_parameters_HPPC.mat','param')
%--------------------------------------------------------------------------
% Option 2: If you don't have identified parameters, then load  parameters from literature [Chen 2020]
%--------------------------------------------------------------------------
% param = Parameters_LG_INR21700_M50;

%--------------------------------------------------------------------------
% Initial SOC [-]
%--------------------------------------------------------------------------
SOC_init = 1;

%--------------------------------------------------------------------------
% Load HPPC profile input current vector [A] and simulation time [s] 
%--------------------------------------------------------------------------
% time_vec: Should be a vector where each entry defines the time period for each CC segment [s] (Px1)
% curr_vec: Should be a vector where each entry defines the current for that CC segment [A] (Px1) (negative current: discharging)
% -> where P is the number of CC segments for your HPPC experiment
%--------------------------------------------------------------------------
% For example, since HPPC starts with a 1C discharge for 360 s followed by
% a rest period for 3600 s:
%   -> time_vec(1) = 360
%   -> curr_vec(1) = -4.8504
%   -> time_vec(2) = 3600
%   -> curr_vec(2) = 0
%--------------------------------------------------------------------------
% Load file that contains curr_vec and time_vec
load('HPPC_W8_Diag1.mat') % Data extrated cell W8 at Diag.#1 from [2] 

%--------------------------------------------------------------------------
% Current mode (do not change)
%--------------------------------------------------------------------------
% param.CurrentMode = 1     -> constant current
% param.CurrentMode = 2     -> variable current, e.g., UDDS
%--------------------------------------------------------------------------
% NOTE: During HPPC, we are using concatenated CC simulations (much faster than variable current)
%--------------------------------------------------------------------------
param.CurrentMode = 1; 

%--------------------------------------------------------------------------
% Print simulation time 
%--------------------------------------------------------------------------
% param.sim_print = 0     -> supress simulation time
% param.sim_print = 1     -> print simulation time
%--------------------------------------------------------------------------
param.sim_print = 1;

%% HPPC simulation

% Define current density vector for each CC segment
currentDensity_vec = curr_vec/param.Acell;

% Number of CC segments (equal to P)
curr_segments = length(curr_vec); 

% Initial simulation time [s]
t0 = 0;
% Final simulation time vector [s]
tf_vec = time_vec; 

% HPPC simulation
[t_sim,V_sim,SOCp_sim,SOCn_sim,curr_dens_sim,...
    delta_phie_sim,delta_eta_sim,delta_OCV_sim,~]=HPPC_sim(param,SOC_init,t0,tf_vec,currentDensity_vec,curr_segments);

%% Visualizing results 
close all

% Voltage-current plot
figure; hold on; box on;
yyaxis left
plot(t_sim./3600,V_sim,'linewidth',2)
xlabel('Time [h]'); ylabel('Voltage, $V$ [V]')
ylim([param.lowerCutoffVoltage param.upperCutoffVoltage])
yyaxis right
plot(t_sim./3600,curr_dens_sim*param.Acell,'linewidth',2)
xlabel('Time [h]'); ylabel('Current, $I$ [A]');xlim([0 inf])
title('HPPC Simulation')
set(gca,'Fontsize',32);
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','latex')
set(gcf,'color','w');
set(gcf, 'Units','inches','Position', [1 1 13 8]);
legend off

% SOC plot
figure; hold on; box on;
plot(t_sim./3600,SOCp_sim*100,'linewidth',2)
plot(t_sim./3600,SOCn_sim*100,'linewidth',2)
xlabel('Time [h]'); ylabel('State-of-charge, $SOC$ [\%]');legend('$SOC_p$','$SOC_n$');xlim([0 inf])
title('HPPC Simulation')
set(gca,'Fontsize',32);
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','latex')
set(gcf,'color','w');
set(gcf, 'Units','inches','Position', [1 1 13 8]);

% Open-circuit potential, overpotential, ohmic drop difference between x=0 (CC|positive) and x=L (negative|CC)
%------------------------------------------------------------------------
% Voltage equation in DFN model is V = phis_p(x=0) - phis_n(x=L) + Rc*I (negative current is discharging),
% where:
    % x = 0 corresponds to interface between cathode current collector and cathode (CC|positive) 
    % x = L corresponds to interaface between anode current collector and anode (negative|CC)
% Since phis_p(x=0) = eta_p(x=0) + phie(x=0) + Up(x=0) and phis_n(x=L) = eta_n(x=L) + phie(x=L) + Un(x=L),
% V = eta_p(x=0) + phie(x=0) + Up(x=0) - (eta_p(x=L) + phie(x=L) + Un(x=L)) + Rc*I
% V = delta_eta + delta_phie + delta_OCV + Rc*I
% where:
    % delta_eta = eta_p(x=0) - eta_n(x=L)
    % delta_phie = phie(x=0) - phie(x=L)
    % delta_OCV = Up(x=0) - Un(x=L)
%------------------------------------------------------------------------
figure; hold on
yyaxis left
plot(t_sim/3600,V_sim,'Linewidth',2,'DisplayName','$V=\Delta OCP+\Delta \phi_e+\Delta \eta + IR_c$')
plot(t_sim/3600,delta_OCV_sim,'Linewidth',2,'DisplayName','$\Delta OCP = U_p(x=0)-U_n(x=L)$')
ylabel('Voltage [V]'); xlabel('Time [h]')
xlabel('Time [s]'); ylabel('[V]')
yyaxis right
plot(t_sim/3600,delta_phie_sim,'Linewidth',2,'DisplayName','$\Delta \phi_e = \phi_e(x=0)-\phi_e(x=L)$')
plot(t_sim/3600,delta_eta_sim,'Linewidth',2,'DisplayName','$\Delta \eta = \eta_p(x=0)-\eta_n(x=L)$')
plot(t_sim/3600,curr_dens_sim*param.Rc*param.Acell,'Linewidth',2,'DisplayName','$IR_c$')
xlabel('Time [s]'); ylabel('[V]');xlim([0 inf]);ylim([-0.2 0.2])
title('HPPC Simulation')
set(gca,'Fontsize',24);legend('location','southeast')
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','latex')
set(gcf,'color','w');
set(gcf, 'Units','inches','Position', [1 1 13 9]);
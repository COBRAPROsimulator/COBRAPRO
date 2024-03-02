% -------------------------------------------------------------------------
% DFN_LSA_CC.m
% -------------------------------------------------------------------------
%   DFN_LSA_CC performs local sensitivity analysis (LSA) by perturbing 
%   each parameter one-at-a-time. This script can be used to determine sensitive
%   parameters for a given CC usage senario. Once the sensitive parameters are
%   determined, they can be identified using CC experimental data using "DFN_pso_CC.m" script.

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
% Load nominal/reference parameters
%--------------------------------------------------------------------------
param = Parameters_LG_INR21700_M50;

%--------------------------------------------------------------------------
% Enter names of parameters to conduct LSA (make sure names match the
% parameter names in "param" structure containing nominal parameters)
%--------------------------------------------------------------------------
param_LSA_CC = {'ep' 'en' 'Rpp' 'Rpn'};

% Type in latex version of parameter names for nice plots
theta_names = {'$\varepsilon_p$' '$\varepsilon_n$' '$R_p$' '$R_n$'};

%--------------------------------------------------------------------------
% Perturbation coefficient for LSA [-]
%--------------------------------------------------------------------------
pct=0.05;

%--------------------------------------------------------------------------
% Define input current density [A.m-2]
%--------------------------------------------------------------------------
% Method 1: Define C-rate (negative C-rate means discharging)
%   -> Crate = -1; 
%   -> param.current = Crate*param.Q_nom;
% Method 2: Define current directly (negative current means discharging)
%   -> param.current = -4.85; 
%--------------------------------------------------------------------------
Crate = -0.05; % [A/Ah]
param.current = Crate*param.Q_nom; % [A]

% Define current density based on input current 
param.currentDensity = param.current/param.Acell; % [A.m-2]

%--------------------------------------------------------------------------
% Initial SOC [-]
%--------------------------------------------------------------------------
param.SOC_init = 1;

%--------------------------------------------------------------------------
% Simulation time 
%--------------------------------------------------------------------------
% NOTE: Simulation will exit automatically when the lower or upper cut-off
% voltages are reached. To simulate till the cut-off voltage condition,
% keep tf large enough, e.g., for C/20 the tf can be greater than 72000s.
%--------------------------------------------------------------------------
% Initial simulation time [s]
param.t0 = 0;
% Final simulation time [s]
param.tf = 1e6;

%--------------------------------------------------------------------------
% Current mode (do not change)
%--------------------------------------------------------------------------
% param.CurrentMode = 1     -> constant current
% param.CurrentMode = 2     -> variable current, e.g., UDDS
%--------------------------------------------------------------------------
param.CurrentMode = 1;

%% Extract nominal parameters

theta_nom = zeros(1,length(theta_names));
for i = 1:length(theta_names)
    theta_nom(i) = param.(param_LSA_CC{i});
end

%% Simulate model for each perturbed parameter 

[t_nom,V_nom,SOCp_nom,SOCn_nom,t_l_all,V_l_all,SOCp_l_all,SOCn_l_all,...
    t_u_all,V_u_all,SOCp_u_all,SOCn_u_all,...
    S_norm_V, S_norm_SOCp, S_norm_SOCn, S_norm_V_SOC]...
    = DFN_LSA_sim_CC(param,param_LSA_CC,theta_nom,pct);

%% Voltage, SOCp, SOCn plots as a function of time

for i = 1:length(theta_names)
    figure; hold on;set(gcf,'color','w')
    set(gcf, 'Units','inches','Position', [1 1 20 8]);
    subplot(1,3,1)
    plot(t_nom,V_nom,'linewidth',2,'DisplayName',['$\theta_{nom}$=' num2str(theta_nom(i))])
    hold on
    plot(t_l_all{i},V_l_all{i},'linewidth',2,'DisplayName',['$\theta_{lower}$=' num2str(theta_nom(i)*(1-pct))])
    plot(t_u_all{i},V_u_all{i},'linewidth',2,'DisplayName',['$\theta_{upper}$=' num2str(theta_nom(i)*(1+pct))])
    title(['Voltage ' theta_names{i} ', $\overline{S_{V}}$:' num2str(vecnorm(S_norm_V{i}))])
    ylim([-inf inf])
    set(gca,'Fontsize',32);legend('location','best')
    set(findall(gcf,'-property','interpreter'),'interpreter','latex')
    set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','latex')
    subplot(1,3,2) 
    plot(t_nom,SOCp_nom,'linewidth',2,'DisplayName',['$\theta_{nom}$=' num2str(theta_nom(i))])
    hold on
    plot(t_l_all{i},SOCp_l_all{i},'linewidth',2,'DisplayName',['$\theta_{lower}$=' num2str(theta_nom(i)*(1-pct))])
    plot(t_u_all{i},SOCp_u_all{i},'linewidth',2,'DisplayName',['$\theta_{upper}$=' num2str(theta_nom(i)*(1+pct))])
    title(['$SOC_p$ ' theta_names{i} ', $\overline{S_{SOC_p}}$:' num2str(vecnorm(S_norm_SOCp{i}))])
    set(gca,'Fontsize',32);legend('location','best')
    set(findall(gcf,'-property','interpreter'),'interpreter','latex')
    set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','latex')
    subplot(1,3,3) 
    plot(t_nom,SOCn_nom,'linewidth',2,'DisplayName',['$\theta_{nom}$=' num2str(theta_nom(i))])
    hold on
    plot(t_l_all{i},SOCn_l_all{i},'linewidth',2,'DisplayName',['$\theta_{lower}$=' num2str(theta_nom(i)*(1-pct))])
    plot(t_u_all{i},SOCn_u_all{i},'linewidth',2,'DisplayName',['$\theta_{upper}$=' num2str(theta_nom(i)*(1+pct))])
    title(['$SOC_n$ ' theta_names{i} ', $\overline{S_{SOC_n}}$:' num2str(vecnorm(S_norm_SOCn{i}))])
    set(gca,'Fontsize',32);legend('location','best')
    set(findall(gcf,'-property','interpreter'),'interpreter','latex')
    set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','latex')
end

%% Sensitivity plots as a function of time

figure; hold on
for i = 1:length(theta_names)
    plot(t_l_all{i}(1:length(S_norm_V{i}))/3600,S_norm_V{i},'LineWidth',2,'DisplayName',theta_names{i})
end
xlabel('Time [h]');ylabel('$S_V$ [-]');xlim([0 inf])
title('Normalized Voltage Sensitivity Vector as a function of time')
set(gca,'Fontsize',32);legend('location','best');set(gcf,'color','w')
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','latex')
set(gcf, 'Units','inches','Position', [1 1 15 8]);

figure; hold on
for i = 1:length(theta_names)
    subplot(2,1,1); hold on
    plot(t_l_all{i}(1:length(S_norm_V{i}))/3600,S_norm_SOCp{i},'LineWidth',2,'DisplayName',theta_names{i})
    subplot(2,1,2); hold on
    plot(t_l_all{i}(1:length(S_norm_V{i}))/3600,S_norm_SOCn{i},'LineWidth',2,'DisplayName',theta_names{i})
end

subplot(2,1,1)
xlabel('Time [h]');ylabel('$S_{SOC_p}$ [-]');xlim([0 inf])
title('Normalized $SOC_p$ Sensitivity Vector as a function of time')
set(gca,'Fontsize',32);legend('location','best');set(gcf,'color','w')
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','latex')
set(gcf, 'Units','inches','Position', [1 1 15 8]);

subplot(2,1,2)
xlabel('Time [h]');ylabel('$S_{SOC_n}$ [-]');xlim([0 inf])
title('Normalized $SOC_n$ Sensitivity Vector as a function of time')
set(gca,'Fontsize',32);legend('location','best');set(gcf,'color','w')
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','latex')
set(gcf, 'Units','inches','Position', [1 1 15 8]);

%% Compute 2-norm of each column, representing the average sensitivity of each parameter

S_2norm_V_all = zeros(1,length(theta_names));
S_2norm_V_SOC_all = zeros(1,length(theta_names));
length_V_SOC_vec = zeros(1,length(theta_names));

for i = 1:length(theta_names)
    % Compute 2-norm of voltage sensitivity vector
    S_2norm_V_all(i) = vecnorm(S_norm_V{i});
    % Compute 2-norm of voltage + SOCp + SOCn sensitivity vector
    S_2norm_V_SOC_all(i) = vecnorm(S_norm_V_SOC{i});
    length_V_SOC_vec(i) = length(S_norm_V_SOC{i});
end

[S_2norm_V_SOC_all_sort,ind] = sort(S_2norm_V_SOC_all,'descend');
S_2norm_V_all_sort = S_2norm_V_all(ind);

%% Normalized sensitivity dot plots for each parameter 

figure; hold on; grid on
plot(S_2norm_V_SOC_all_sort,'o','MarkerSize',12,'linewidth',2,'DisplayName','$|S_{V,SOC}|$')
plot(S_2norm_V_all_sort,'^','MarkerSize',12,'linewidth',2,'DisplayName','$|S_{V}|$')
xticks(1:length(theta_names));xticklabels(theta_names(ind));
ylabel('$|S|$ [-]')
set(gca, 'YScale', 'log')
title('2-norm of Normalized Sensitivity Vector');legend('location','best')
set(gcf,'color','w')
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','latex')
set(gcf, 'Units','inches','Position', [1 1 15 8]);set(gca,'Fontsize',32);
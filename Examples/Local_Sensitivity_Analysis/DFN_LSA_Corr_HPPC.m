% -------------------------------------------------------------------------
% DFN_LSA_Corr_HPPC.m
% -------------------------------------------------------------------------
% DFN_LSA_Corr_HPPC performs local sensitivity analysis (LSA) for the HPPC profile
% by perturbing each parameter one at a time. Also, correlation analysis is conducted
% to quantify parameter correlation. Script will output identifiable
% parameters using just LSA and using LSA correlation analysis. 
% The sensitive/uncorrelated parameters can be identified using HPPC
% experimental data in the "DFN_pso_HPPC.m" script.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COBRAPRO: Co-simulation Battery Modeling for Accelerated Parameter Optimization

% Copyright (c) 2024 CO-simulation BatteRy modeling for Accelerated PaRameter Optimization (COBRAPRO)
% COBRAPRO is freely distributed software under the MIT License 
% v1.0.0.: Released 02/2024 

% Written by Sara Ha (sungyeon.sara.ha@stanford.edu)
% PI: Prof. Simona Onori (sonori@stanford.edu)
% Stanford Energy Control Group, Stanford University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%% User Input  
clc;clear;close all

%--------------------------------------------------------------------------
% Option 1: Enter identified parameters
%--------------------------------------------------------------------------
load('identified_parameters_0_05C.mat','param')
param.upperCutoffVoltage = 4.3; 
%--------------------------------------------------------------------------
% Option 2: If you don't have identified parameters, then load parameters from literature [Chen 2020]
%--------------------------------------------------------------------------
% param = Parameters_LG_INR21700_M50;

%--------------------------------------------------------------------------
% Option 1: Enter names of parameters to conduct LSA (make sure names match the
% parameter names in "param" structure containing nominal parameters)
%--------------------------------------------------------------------------
% param_LSA_HPPC = {'Dsp' 'Dsn' 't1_constant' 'kp' 'kn' 'c0' 'De' 'Kappa'};
%--------------------------------------------------------------------------
% Option 2: Load identifiable parameters from identifiability analysis
% conducted in "Examples/Local_Sensitivity_Analysis/DFN_LSA_Corr_HPPC.m"
%--------------------------------------------------------------------------
load('DFN_identification_results/HPPC_identifiable_params.mat','LSA_identifiable','corr_identifiable')
% 1. Load results from LSA analysis 
param_LSA_HPPC = LSA_identifiable;
% 2. Load results from LSA + correlation analysis 
% param_LSA_HPPC = corr_identifiable;

% Type in latex version of param_LSA_HPPC for nice plots
theta_names = {'$D_{s,p}$' '$k_n$' '$c_0$' '$D_{s,n}$'};

%--------------------------------------------------------------------------
% Perturbation coefficient for LSA [-]
%--------------------------------------------------------------------------
pct=0.05;

%--------------------------------------------------------------------------
% Define sensitivity index threshold
%--------------------------------------------------------------------------
% Parameters with sensitivity index higher than beta_LSA will be considered
% the "identifiable" parameters for the HPPC profile
%--------------------------------------------------------------------------
beta_LSA = 0.5;

%--------------------------------------------------------------------------
% Define correlation coefficient threshold
%--------------------------------------------------------------------------
% Correlation coefficient threshold used to determine the uncorrelated parameter  
% vector, which will be considered the "identifiable" parameters for the HPPC profile
%--------------------------------------------------------------------------
beta_corr = 0.8;

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
% Initial SOC [-]
%--------------------------------------------------------------------------
param.SOC_init = 1;

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
    theta_nom(i) = param.(param_LSA_HPPC{i});
end

%% Processing simulation settings

% Number of CC segments
curr_segments = length(curr_vec); 

% Initial simulation time [s]
param.t0 = 0;
% Final simulation time vector [s]
param.tf_vec = time_vec;

% Add simulation settings to param for PSO
param.curr_vec = curr_vec;
param.curr_dens_vec = curr_vec/param.Acell;
param.curr_segments = curr_segments;

%% Simulate model for each perturbed parameter 

[t_nom,V_nom,SOCp_nom,SOCn_nom,t_l_all,V_l_all,SOCp_l_all,SOCn_l_all,...
    t_u_all,V_u_all,SOCp_u_all,SOCn_u_all,...
    S_norm_V, S_norm_SOCp, S_norm_SOCn, S_norm_V_SOC]...
    = DFN_LSA_sim_HPPC(param,param_LSA_HPPC,theta_nom,pct);

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
colorful = brewermap([],'Paired');
figure; hold on
for i = 1:length(theta_names)
    plot(t_l_all{i}(1:length(S_norm_V{i}))/3600,S_norm_V{i},'color',colorful(i,:),'LineWidth',2,'DisplayName',theta_names{i})
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
    plot(t_l_all{i}(1:length(S_norm_V{i}))/3600,S_norm_SOCp{i},'LineWidth',2,'color',colorful(i,:),'DisplayName',theta_names{i})
    subplot(2,1,2); hold on
    plot(t_l_all{i}(1:length(S_norm_V{i}))/3600,S_norm_SOCn{i},'LineWidth',2,'color',colorful(i,:),'DisplayName',theta_names{i})
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
plot(S_2norm_V_SOC_all_sort,'o','MarkerSize',12,'linewidth',2,'DisplayName','$S^{V,SOC}$')
plot(S_2norm_V_all_sort,'^','MarkerSize',12,'linewidth',2,'DisplayName','$S^{V}$')
xticks(1:length(theta_names));xticklabels(theta_names(ind));
ylabel('Sensitivity [-]')
set(gca, 'YScale', 'log'); 
title('HPPC Sensitivity');%title('2-norm of Normalized Sensitivity Vector')
set(gcf,'color','w');legend('location','best')
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','latex')
set(gcf, 'Units','inches','Position', [1 1 14 8]);set(gca,'Fontsize',32);

%% Correlation matrix calculation

 min_length_V_SOC_vec = min(length_V_SOC_vec);

 S_2norm_V_SOC_matrix = zeros(min_length_V_SOC_vec,length(theta_names));
 for i=1:length(theta_names)
     S_2norm_V_SOC_matrix(:,i) = S_norm_V_SOC{i}(1:min_length_V_SOC_vec);
 end
corr_V_SOC_matrix = corrcoef(S_2norm_V_SOC_matrix);

%% Correlation matrix plot

plot_matrix=corr_V_SOC_matrix;
plot_matrix=abs(plot_matrix);
figure;
h = heatmap(plot_matrix, 'colormap',flipud(parula), 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ");
h.NodeChildren(3).XAxis.TickLabelInterpreter = 'latex';
h.NodeChildren(3).YAxis.TickLabelInterpreter = 'latex';
h.XDisplayLabels = (theta_names);
h.YDisplayLabels = (theta_names);
set(gcf,'color','w');
title('$S_{V,SOC}$ Correlation Matrix')
set(findall(gcf,'-property','FontSize'),'FontSize',20);
set(findall(gcf,'-property','FontNameinterpreter'),'latex')
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','latex')

%% Output and save the identifiable parameters

% When using just LSA to determine identifiable parameters:
param_LSA_HPPC_sorted = param_LSA_HPPC(ind);
LSA_identifiable = param_LSA_HPPC_sorted(S_2norm_V_SOC_all_sort>beta_LSA);

% When using also correlation anlaysis to determine identifiable parameters:
corr_identifiable = unCorr_parameters(param_LSA_HPPC, S_2norm_V_SOC_all, corr_V_SOC_matrix, beta_corr);

% Print results
fprintf('------------------------------------------------\n')
fprintf('Identifiable parameters from LSA analysis:\n')
fprintf(['(sensitivity threshold = ' num2str(beta_LSA) ')\n'])
fprintf('------------------------------------------------\n')
fprintf('%s ',LSA_identifiable{:})
fprintf('\n\n------------------------------------------------\n')
fprintf('Identifiable parameters from LSA and correlation analysis:\n')
fprintf(['(correlation coefficient threshold = ' num2str(beta_corr) ')\n'])
fprintf('------------------------------------------------\n')
fprintf('%s ',corr_identifiable{:})
fprintf('\n')

% Save results
save('DFN_identification_results/HPPC_identifiable_params.mat','LSA_identifiable','corr_identifiable','beta_LSA','beta_corr')
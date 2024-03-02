% -------------------------------------------------------------------------
% DFN_pso_0_05C_validation.m
% -------------------------------------------------------------------------
% DFN_pso_0_05C_validation shows the parameter identification results using
% C/20 experimental data.
%
% Parameters identified during C/20 discharge: {theta100_p,theta0_p,theta100_n,theta0_n}

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
% Enter your identified parameters 
%--------------------------------------------------------------------------
load('identified_parameters_0_05C.mat')

%% Run simulation

current = mean(I_exp);
currentDensity = current/param.Acell;

x_init = [];
results = runModel(x_init,param,t0,tf,SOC_init,currentDensity);

%% Print identified values

fprintf('Displaying identified values...\n')
for i = 1:length(param_CC)
    fprintf('------------------------\n')
    fprintf([param_CC{i} ':\n'])
    fprintf(['Identified value: ' num2str(x_opt_CC(i)) '\n'])
    fprintf([num2str(lb(i)) '(lower) | ' num2str(init_pos(i)) '(initial) | ' num2str(ub(i)) '(upper)\n'])
end

if param.csmax_constraint == 1
    fprintf('\nStoichiometric identification: maximum solid concentration constraint was enabled\n')
    param.ctp = 3600*param.Q_nom/(param.ep_s*param.F*param.lp*param.Acell*(param.theta0_p-param.theta100_p));
    param.ctn = 3600*param.Q_nom/(param.en_s*param.F*param.ln*param.Acell*(param.theta100_n-param.theta0_n));
    fprintf(['csp_max = ' num2str(param.ctp) '\n'])
    fprintf(['csn_max = ' num2str(param.ctn) '\n'])
end

%% Visualization: plot identified results against experimental data

% Simulated data
t_sim = results.t;
V_sim = results.voltage;
SOCp_sim = results.soc_p_bulk;
SOCn_sim = results.soc_n_bulk;

figure; hold on; box on
myexp_color = [0.5138 0.7332 0.8587];

subplot(2,1,1);hold on
plot(t_exp/3600,SOCp_exp,'color',myexp_color,'Linewidth',3,'DisplayName','Experiment')
plot(t_sim/3600,SOCp_sim,'k--','Linewidth',3,'DisplayName','Simulation')
xlabel('Time [h]');ylabel('$SOC_p$ [-]');xlim([0 inf])
set(gca,'Fontsize',32);
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','latex')

subplot(2,1,2);hold on
plot(t_exp/3600,SOCn_exp,'color',myexp_color,'Linewidth',3,'DisplayName','Experiment')
plot(t_sim/3600,SOCn_sim,'k--','Linewidth',3,'DisplayName','Simulation')
xlabel('Time [h]');ylabel('$SOC_n$ [-]');xlim([0 inf])
set(gca,'Fontsize',32);
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','latex')
set(gcf,'color','w');set(gcf, 'Units','inches','Position', [1 1 13 8]);

figure; hold on; box on
plot(t_exp/3600,V_exp,'color',myexp_color,'Linewidth',4,'DisplayName','Experiment')
plot(t_sim/3600,V_sim,'k--','Linewidth',4,'DisplayName','Simulation')
xlabel('Time [h]'); ylabel('Voltage [V]')
ylim([param.lowerCutoffVoltage param.upperCutoffVoltage]);xlim([0 inf])
set(gca,'Fontsize',32);
legend('Experiment','Simulation','Location','NorthEast')
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','latex')
set(gcf,'color','w');set(gcf, 'Units','inches','Position', [1 1 13 8]);

% Interpolate sim results to match the timestep of voltage experiment
V_sim_obj = interp1(t_sim, V_sim, t_exp);
SOCp_obj = interp1(t_sim, SOCp_sim, t_exp);
SOCn_obj = interp1(t_sim, SOCn_sim, t_exp);
% Remove NaN from interpolated data
V_sim_obj=V_sim_obj(~isnan(V_sim_obj));
SOCp_obj=SOCp_obj(~isnan(SOCp_obj));
SOCn_obj=SOCn_obj(~isnan(SOCn_obj));
% Use shorter length to compare
len1=length(V_sim_obj);
len2=length(V_exp);
len=min(len1,len2);

% Objective function
J1 = rms((V_exp(1:len)-V_sim_obj(1:len))./V_sim_obj(1:len));
J1_mV = rms(V_exp(1:len)-V_sim_obj(1:len));
J2 = rms(SOCp_exp(1:len)-SOCp_obj(1:len));
J3 = rms(SOCn_exp(1:len)-SOCn_obj(1:len));
J_tot=J1+J2+J3;

% Print objective function values
fprintf('\nDisplaying objective function values...\n')
fprintf('------------------------\n')
fprintf(['J_V =' num2str(J1) ' [-]\n'])
fprintf(['J_V =' num2str(J1_mV*1000) ' [mV]\n'])
fprintf(['J_SOCp =' num2str(J2*100) ' [%%]\n'])
fprintf(['J_SOCn =' num2str(J3*100) ' [%%]\n'])
fprintf(['J_tot =' num2str(J1+J2+J3) ' [-]\n'])
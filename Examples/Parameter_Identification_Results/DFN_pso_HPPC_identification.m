% -------------------------------------------------------------------------
% DFN_pso_HPPC_validation.m
% -------------------------------------------------------------------------
% DFN_pso_HPPC_validation shows the parameter identification results using
% HPPC experimental data. Prior to HPPC identification, C/20 capacity test
% data was used to identify the stoichiometric parameters.
%
% Parameters identified during C/20 discharge: {theta100_p,theta0_p,theta100_n,theta0_n}
% Parameters identified during HPPC: {kn, Dsp, c0, kp, Dsn, Kappa, De, t1_constant}

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
load('identified_parameters_HPPC.mat')

%--------------------------------------------------------------------------
% Print simulation time 
%--------------------------------------------------------------------------
% param.sim_print = 0     -> supress simulation time
% param.sim_print = 1     -> print simulation time
%--------------------------------------------------------------------------
param.sim_print = 1;

%% HPPC simulation

% HPPC simulation
[t_sim,V_sim,SOCp_sim,SOCn_sim,curr_dens_sim,...
    delta_phie_sim,delta_eta_sim,delta_OCV_sim,~]=HPPC_sim(param,SOC_init,t0,tf_vec,currentDensity_vec,curr_segments);

%% Print identified values
fprintf('\nDisplaying identified values...\n')
for i = 1:length(param_HPPC)
    fprintf('------------------------\n')
    fprintf([param_HPPC{i} ':\n'])
    fprintf(['Identified value: ' num2str(x_opt_HPPC(i)) '\n'])
    fprintf([num2str(lb(i)) '(lower) | ' num2str(init_pos(i)) '(initial) | ' num2str(ub(i)) '(upper)\n'])
end

%% Visualizing results 

figure; hold on; box on; set(gcf,'color','w');
myexp_color = [0.5138 0.7332 0.8587];

subplot(2,1,1)
plot(t_exp/3600,SOCp_exp,'color',myexp_color,'Linewidth',3,'DisplayName','Experiment')
hold on
plot(t_sim/3600,SOCp_sim,'k--','Linewidth',3,'DisplayName','Simulation')
xlabel('Time [h]');ylabel('$SOC_p$ [-]');xlim([0 inf])
set(gca,'Fontsize',32);
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','latex')

subplot(2,1,2)
plot(t_exp/3600,SOCn_exp,'color',myexp_color,'Linewidth',3,'DisplayName','Experiment')
hold on
plot(t_sim/3600,SOCn_sim,'k--','Linewidth',3,'DisplayName','Simulation')
xlabel('Time [h]');ylabel('$SOC_n$ [-]');xlim([0 inf])
set(gca,'Fontsize',32);
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','latex')
set(gcf, 'Units','inches','Position', [1 1 13 8]);

figure; hold on; box on
plot(t_exp/3600,V_exp,'color',myexp_color,'Linewidth',4,'DisplayName','Experiment')
plot(t_sim/3600,V_sim,'k--','Linewidth',4,'DisplayName','Simulation')
xlabel('Time [h]'); ylabel('Voltage [V]')
ylim([param.lowerCutoffVoltage param.upperCutoffVoltage]);xlim([0 inf])
legend('Experiment','Simulation','Location','NorthEast')
set(gca,'Fontsize',32);
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','latex')
set(gcf,'color','w');
set(gcf, 'Units','inches','Position', [1 1 13 8]);

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
fprintf(['J_SOCp =' num2str(J2*100)  ' [%%]\n'])
fprintf(['J_SOCn =' num2str(J3*100)  ' [%%]\n'])
fprintf(['J_tot =' num2str(J1+J2+J3)  ' [-]\n'])
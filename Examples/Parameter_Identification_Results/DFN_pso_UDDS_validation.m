% -------------------------------------------------------------------------
% DFN_pso_UDDS_validation.m
% -------------------------------------------------------------------------
% DFN_pso_UDDS_validation shows UDDS validation identified parameters from
% C/20 discharge and HPPC data.
%
% Parameters identified during C/20 discharge: {theta100_p,theta0_p,theta100_n,theta0_n}
% Parameters identified during HPPC: {kp, kn, Dsp, Dsn, Kappa, De, t1_constant, c0}

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
clear;clc;close all

%--------------------------------------------------------------------------
% Load identification results 
%--------------------------------------------------------------------------
load('identified_parameters_HPPC.mat')

%--------------------------------------------------------------------------
% Load Experimental Data 
%--------------------------------------------------------------------------
%   t: Should be a vector consisting of your time experiment data      [s] (Mx1)
%   I: Should be a vector consisting of your current experiment data   [A] (Mx1) (negative current: discharging)
%   V: Should be a vector consisting of your volatge experiemntal data [V] (Mx1)
%   -> where M is the total number of data points in your experiment
%--------------------------------------------------------------------------
% HPPC test conducted on LG INR21700 M50T cells
load('UDDS_W8_cyc1.mat')

t = t_data';
I = I_data';
V = V_data';

%--------------------------------------------------------------------------
% Initial SOC [-]
%--------------------------------------------------------------------------
SOC_init = 0.7816; 

%--------------------------------------------------------------------------
% Simulation time 
%--------------------------------------------------------------------------
% Initial simulation time [s]
t0 = 0;
% Final simulation time [s]
tf = t_data(end);

%--------------------------------------------------------------------------
% Choose the desired sampling time [s] to resample your experiment data 
%--------------------------------------------------------------------------
deltaT_exp = 1; % [s]

%--------------------------------------------------------------------------
% Current mode (do not change)
%--------------------------------------------------------------------------
% param.CurrentMode = 1     -> constant current
% param.CurrentMode = 2     -> variable current, e.g., UDDS
%--------------------------------------------------------------------------
param.CurrentMode = 2;

%--------------------------------------------------------------------------
% IDA solver type (do not change)
%--------------------------------------------------------------------------
% param.IDA_type = 1 -> OneStep (IDA solver determines efficient time step--usually faster than Normal type)
% param.IDA_type = 2 -> Normal  (user defined time step)
%--------------------------------------------------------------------------
% For variable current, we found that it is more efficient to use Normal type
param.IDA_type = 2;
param.deltaT = 1; % [s] (only for Normal type)

%--------------------------------------------------------------------------
% Print simulation time
%--------------------------------------------------------------------------
% param.sim_print = 0     -> supress simulation time
% param.sim_print = 1     -> print simulation time
%--------------------------------------------------------------------------
param.sim_print = 1;

%% Processing experimental data 

% Interpolate experimental data according to user defined sampling time
t_exp = (0:deltaT_exp:t(end))'; 
I_exp = interp1(t,I,t_exp);     % Negative current: discharging
V_exp = interp1(t,V,t_exp);

% Calculate SOC using coulomb counting
SOC_exp = SOC_init - cumtrapz(t_exp/3600,-I_exp)/param.Q_nom;

% Ensure that SOC is not below 0% or above 100% 
SOC_exp(SOC_exp<0) = 0;
SOC_exp(SOC_exp>1) = 1;

% Define SOCp and SOCn
SOCp_exp = SOC_exp;
SOCn_exp = SOC_exp;

%% UDDS simulation

% For variable current mode, currentDensity is a dummy variable and not used in simulation. Set to zero.
currentDensity = 0;
% Pass the UDDS time and current vectors to param 
param.t_data = t_exp;
param.I_data = I_exp;

% For variable current, must define the getCurrentDensity.m function as param.currentDensityFunction
% -> getCurrentDensity.m calculates time varying current density based on param.t_data and param.I_data
param.currentDensityFunction = @getCurrentDensity;

x_init = [];
%--------------------------------------------------------------------------
% Option 1: Run UDDS simulation (takes a few minutes to complete since variable current and UDDS simulation time is long)
%--------------------------------------------------------------------------
% [output,param] = runModel(x_init,param,t0,tf,SOC_init,currentDensity);
%--------------------------------------------------------------------------
% Option 2: To save time, load pregenerated UDDS simulation results (will yield same results as Option 1)
%--------------------------------------------------------------------------
load('UDDS_validation_results.mat')

%% Visualizing results 

t_sim = output.t;
V_sim = output.voltage;
SOCp_sim = output.soc_p_bulk;
SOCn_sim = output.soc_n_bulk;

myexp_color = [0.5138 0.7332 0.8587];

figure; hold on; box on; set(gcf,'color','w');
set(gcf, 'Units','inches','Position', [1 1 13 8]);
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

figure; hold on; box on;
plot(t_exp,V_exp,'color',myexp_color,'linewidth',2,'DisplayName','Experimental')
plot(output.t,output.voltage,'k','linewidth',2,'DisplayName','Simulation')
xlabel('Time [h]'); ylabel('Voltage [V]')
ylim([3.2 4.1]);xlim([0 inf])
legend('location','best')
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
J1 = rms((V_exp(1:len)-V_sim_obj(1:len))./V_exp(1:len));
J1_mV = rms(V_exp(1:len)-V_sim_obj(1:len));
J2 = rms(SOC_exp(1:len)-SOCp_obj(1:len));
J3 = rms(SOC_exp(1:len)-SOCn_obj(1:len));
J_tot=J1+J2+J3;

% Print objective function values
fprintf('\nDisplaying objective function values...\n')
fprintf('------------------------\n')
fprintf(['J_V =' num2str(J1) ' [-]\n'])
fprintf(['J_V =' num2str(J1_mV*1000) ' [mV]\n'])
fprintf(['J_SOCp =' num2str(J2*100)  ' [%%]\n'])
fprintf(['J_SOCn =' num2str(J3*100)  ' [%%]\n'])
fprintf(['J_tot =' num2str(J1+J2+J3)  ' [-]\n'])

% -------------------------------------------------------------------------
% cycle_UDDS.m
% -------------------------------------------------------------------------
% cycle_UDDS is an example code showing how to simulate UDDS driving
% profile using the variable current solver option.

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
% Option 1: Enter identified parameters
%--------------------------------------------------------------------------
load('identified_parameters_HPPC.mat','param')
%--------------------------------------------------------------------------
% Option 2: If you don't have identified parameters, then load parameters from literature [Chen 2020]
%--------------------------------------------------------------------------
% param = Parameters_LG_INR21700_M50;

%--------------------------------------------------------------------------
% Initial SOC [-]
%--------------------------------------------------------------------------
SOC_init = 0.7816; 

%--------------------------------------------------------------------------
% Load UDDS current profile 
%--------------------------------------------------------------------------
%   t: Should be a vector consisting of UDDS time      [s] (Mx1)
%   I: Should be a vector consisting of UDDS current   [A] (Mx1) (negative current: discharging)
%   -> where M is the total number of data points 
%--------------------------------------------------------------------------
% HPPC test conducted on LG INR21700 M50T cells
load('UDDS_W8_cyc1.mat')

t = t_data';
I = I_data';

%--------------------------------------------------------------------------
% Simulation time 
%--------------------------------------------------------------------------
% NOTE: Simulation will exit automatically when the lower or upper cut-off
% voltages are reached. 
%--------------------------------------------------------------------------
% Initial simulation time [s]
t0 = 0;
% Final simulation time [s]
tf = 300;

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
% For variable current, we found that it is typically more efficient to use Normal type
param.IDA_type = 2;
param.deltaT = 1; % [s] (only for Normal type)

%--------------------------------------------------------------------------
% Print simulation time 
%--------------------------------------------------------------------------
% param.sim_print = 0     -> supress simulation time
% param.sim_print = 1     -> print simulation time
%--------------------------------------------------------------------------
param.sim_print = 1;

%% UDDS simulation

% For variable current mode, currentDensity is a dummy variable and not used in simulation. Set to zero.
currentDensity = 0;
% Pass the UDDS time and current vectors to param 
param.t_data = t;
param.I_data = I;

% For variable current, must define the getCurrentDensity.m function as param.currentDensityFunction
% -> getCurrentDensity.m calculates time varying current density based on param.t_data and param.I_data
param.currentDensityFunction = @getCurrentDensity;

% Run UDDS simulation
x_init = [];
[output,param] = runModel(x_init,param,t0,tf,SOC_init,currentDensity);

time = output.t;
voltage = output.voltage;
curr_dens = output.curr_dens;
SOCp = output.soc_p_bulk;
SOCn = output.soc_n_bulk;
ce = output.ce;
csp_surf = output.csp_surf;
csn_surf = output.csn_surf;
phie = output.phie;
etap = output.etap;
etan = output.etan;
Up = output.Up;
Un = output.Un;

%% Visualizing results 
close all

% Voltage-current plot
figure; hold on; box on;
yyaxis left
plot(time,voltage,'linewidth',2)
xlabel('Time [s]'); ylabel('Voltage, $V$ [V]')
ylim([param.lowerCutoffVoltage param.upperCutoffVoltage])
yyaxis right
plot(time,curr_dens*param.Acell,'linewidth',2)
xlabel('Time [s]'); ylabel('Current, $I$ [A]');ylim([-20 20])
title('CC Cycling at 1C')
set(gca,'Fontsize',24);
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','latex')
set(gcf,'color','w');
set(gcf, 'Units','inches','Position', [1 1 13 8]);
legend off

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
phien_ind = param.Np+param.Ns+1:param.Np+param.Ns+param.Nn;
% Calculate delta_phie = phie(x=0) - phie(x=L)
phie_pcc = electrodeBC_hermite_interp(phie(:,1:param.Np),curr_dens,param,'p');
phie_ncc = electrodeBC_hermite_interp(phie(:,phien_ind),curr_dens,param,'n');
delta_phie = phie_pcc - phie_ncc;
% Calculate delta_eta = eta_p(x=0) - eta_n(x=L)
eta_pcc = electrodeBC_hermite_interp(etap,curr_dens,param,'p');
eta_ncc = electrodeBC_hermite_interp(etan,curr_dens,param,'n');
delta_eta = eta_pcc - eta_ncc;
% Calculate delta_OCV = Up(x=0) - Un(x=L)
Upcc = electrodeBC_hermite_interp(Up,curr_dens,param,'p');
Uncc = electrodeBC_hermite_interp(Un,curr_dens,param,'n');
delta_OCV = Upcc - Uncc;

figure; hold on; box on;
yyaxis left
plot(time,voltage,'linewidth',2,'DisplayName','$V=\Delta OCP+\Delta \phi_e+\Delta \eta + IR_c$')
plot(time,delta_OCV,'linewidth',2,'DisplayName','$\Delta OCP = U_p(x=0)-U_n(x=L)$')
xlabel('Time [s]'); ylabel('[V]');ylim([2.5 4.3])
yyaxis right
plot(time,delta_phie,'linewidth',2,'DisplayName','$\Delta \phi_e = \phi_e(x=0)-\phi_e(x=L)$')
plot(time,delta_eta,'linewidth',2,'DisplayName','$\Delta \eta = \eta_p(x=0)-\eta_n(x=L)$')
plot(time,curr_dens*param.Acell*param.Rc,'linewidth',2,'DisplayName','$IR_c$')
xlabel('Time [s]'); ylabel('[V]'); ylim([-1 1])
title('CC Cycling at 1C')
set(gca,'Fontsize',24);legend('location','southeast')
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','latex')
set(gcf,'color','w');
set(gcf, 'Units','inches','Position', [1 1 13 9]);

% State-of-charge in positive and negative electrodes
figure; hold on; box on;
plot(time,SOCp*100,'linewidth',2)
plot(time,SOCn*100,'linewidth',2)
xlabel('Time [s]'); ylabel('State-of-charge, $SOC$ [\%]')
title('CC Cycling at 1C');legend({'$SOC_p$' '$SOC_n$'}) 
set(gca,'Fontsize',24);
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','latex')
set(gcf,'color','w');
set(gcf, 'Units','inches','Position', [1 1 13 8]);

% Electrolyte concentration in positive, separator, negative at t=0,tf/2,tf
figure; hold on; box on;
plot(param.x_intval_psn*1e6,ce(1,:),'k','linewidth',2)
plot(param.x_intval_psn*1e6,ce(round(length(time)/2),:),'r','linewidth',2)
plot(param.x_intval_psn*1e6,ce(length(time),:),'b','linewidth',2)
xline(param.lp*1e6);xline((param.lp+param.ls)*1e6)
xlabel('x-direction [$\mu m$]'); ylabel('Electrolyte concentration, $c_e$ [mol.m-3]');xlim([0 inf])
legend({'t=0 s' ['t=' num2str(round(length(time)/2)) ' s'] ['t=' num2str(length(time)) ' s']}, 'location','southwest')
title('Positive Electrode$|$Separator$|$Negative Electrode')
set(gca,'Fontsize',24);
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','latex')
set(gcf,'color','w');
set(gcf, 'Units','inches','Position', [1 1 13 8]);

% Solid surface concentration in positive and negative electrodes at t=0,tf/2,tf
figure; hold on; box on;
plot(param.x_intval_p*1e6,csp_surf(1,:),'k','linewidth',2)
plot(param.x_intval_p*1e6,csp_surf(round(length(time)/2),:),'r','linewidth',2)
plot(param.x_intval_p*1e6,csp_surf(length(time),:),'b','linewidth',2)
plot(param.x_intval_n*1e6,csn_surf(1,:),'k','linewidth',2)
plot(param.x_intval_n*1e6,csn_surf(round(length(time)/2),:),'r','linewidth',2)
plot(param.x_intval_n*1e6,csn_surf(length(time),:),'b','linewidth',2)
xline(param.lp*1e6);xline((param.lp+param.ls)*1e6)
xlabel('x-direction [$\mu m$]'); ylabel('Solid surface concentration, $c_{s}^{surf}$ [mol.m-3]');xlim([0 inf])
legend({'t=0 s' ['t=' num2str(round(length(time)/2)) ' s'] ['t=' num2str(length(time)) ' s']}, 'location','southwest')
title('Positive Electrode$|$Separator$|$Negative Electrode')
set(gca,'Fontsize',24);
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','latex')
set(gcf,'color','w');
set(gcf, 'Units','inches','Position', [1 1 13 8]);

% Open-circuit potentials
figure; hold on; box on;
plot(param.x_intval_p*1e6,Up(1,:),'k','linewidth',2)
plot(param.x_intval_p*1e6,Up(round(length(time)/2),:),'r','linewidth',2)
plot(param.x_intval_p*1e6,Up(length(time),:),'b','linewidth',2)
plot(param.x_intval_n*1e6,Un(1,:),'k','linewidth',2)
plot(param.x_intval_n*1e6,Un(round(length(time)/2),:),'r','linewidth',2)
plot(param.x_intval_n*1e6,Un(length(time),:),'b','linewidth',2)
xline(param.lp*1e6);xline((param.lp+param.ls)*1e6)
xlabel('x-direction [$\mu m$]'); ylabel('Open-circuit voltage, $OCV$ [V]');xlim([0 inf])
legend({'t=0 s' ['t=' num2str(round(length(time)/2)) ' s'] ['t=' num2str(length(time)) ' s']}, 'location','southwest')
title('Positive Electrode$|$Separator$|$Negative Electrode')
set(gca,'Fontsize',24);
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','latex')
set(gcf,'color','w');
set(gcf, 'Units','inches','Position', [1 1 13 8]);

% Overpotentials
figure; hold on; box on;
plot(param.x_intval_p*1e6,etap(1,:),'k','linewidth',2)
plot(param.x_intval_p*1e6,etap(round(length(time)/2),:),'r','linewidth',2)
plot(param.x_intval_p*1e6,etap(length(time),:),'b','linewidth',2)
plot(param.x_intval_n*1e6,etan(1,:),'k','linewidth',2)
plot(param.x_intval_n*1e6,etan(round(length(time)/2),:),'r','linewidth',2)
plot(param.x_intval_n*1e6,etan(length(time),:),'b','linewidth',2)
xline(param.lp*1e6);xline((param.lp+param.ls)*1e6)
xlabel('x-direction [$\mu m$]'); ylabel('Overpotential, $\eta$ [V]');xlim([0 inf])
legend({'t=0 s' ['t=' num2str(round(length(time)/2)) ' s'] ['t=' num2str(length(time)) ' s']}, 'location','southwest')
title('Positive Electrode$|$Separator$|$Negative Electrode')
set(gca,'Fontsize',24);
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','latex')
set(gcf,'color','w');
set(gcf, 'Units','inches','Position', [1 1 13 8]);
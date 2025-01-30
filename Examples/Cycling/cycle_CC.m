% -------------------------------------------------------------------------
% cycle_CC.m
% -------------------------------------------------------------------------
% cycle_CC is an example code showing how to simulate constant current (CC) cycling.

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
load('identified_parameters_HPPC_0_99Corr.mat','param')
%--------------------------------------------------------------------------
% Option 2: If you don't have identified parameters, then load parameters from literature [Chen 2020]
%--------------------------------------------------------------------------
% param = Parameters_LG_INR21700_M50;

%--------------------------------------------------------------------------
% Number of cycles (one full cycle, e.g. discharging and charging, is equal to two cycles)
%--------------------------------------------------------------------------
cycle_num = 3; 

%--------------------------------------------------------------------------
% Define initial SOC [-] (ranges from 0 to 1)
%--------------------------------------------------------------------------
SOC_init = 1;

%--------------------------------------------------------------------------
% Define input current density [A.m-2]
%--------------------------------------------------------------------------
% Method 1: Define C-rate (negative C-rate means discharging)
%   -> Crate = -1; 
%   -> current = Crate*param.Q_nom;
% Method 2: Define current directly (negative current means discharging)
%   -> current = -4.85; 
%--------------------------------------------------------------------------
Crate = -1; % [A/Ah]
current = Crate*param.Q_nom; % [A]

% Define current density based on input current 
currentDensity = current/param.Acell; % [A.m-2]

%--------------------------------------------------------------------------
% Simulation time [s] for each cycle 
%--------------------------------------------------------------------------
% NOTE: Simulation will exit automatically when the lower or upper cut-off
% voltages are reached. To simulate till the cut-off voltage condition,
% keep tf large enough, e.g., for 1C the tf can be greater than 3600 s.
%--------------------------------------------------------------------------
% Initial simulation time [s]
t0 = 0;
% Final simulation time for each cycle [s]
tf = 4000;

%--------------------------------------------------------------------------
% Current mode (do not change)
%--------------------------------------------------------------------------
% param.CurrentMode = 1     -> constant current
% param.CurrentMode = 2     -> variable current, e.g., UDDS
%--------------------------------------------------------------------------
param.CurrentMode = 1;

%--------------------------------------------------------------------------
% Print simulation time 
%--------------------------------------------------------------------------
% param.sim_print = 0     -> supress simulation time
% param.sim_print = 1     -> print simulation time
%--------------------------------------------------------------------------
param.sim_print = 1;

%% Cycling simulation

% Variables we want to extract from output (for other available interval variables refer to "runModel.m" function)
time = [];
voltage = [];
ce = [];
phie = [];
phisp = [];
phisn = [];
csp_surf = [];
csn_surf = [];
SOCp = [];
SOCn = [];
etap = [];
etan = [];
Up = [];
Un = [];
curr_dens = [];
DAE_solve_time = zeros(1,cycle_num);

% For an empty x_init vector, runModel.m will calculate the initial guess
% based on OCV conditions and the initializer will find the correct
% algebraic initial conditions
x_init = [];

% Start CC cycling
for i = 1:cycle_num
    % Run model for each cyle
    [output,param] = runModel(x_init,param,t0,tf,SOC_init,currentDensity);
    DAE_solve_time(i) = output.DAE_solve_time;
    % Switch current density sign for next cycle
    currentDensity = -currentDensity;
    % Store initial conditions for next cycle
    x_init = output.x_initials;
    % Concatenate output results
    if i == 1
        time(1+size(time,1):size(output.t,1)+size(time,1),1) = output.t;
    else
        time(1+size(time,1):size(output.t,1)+size(time,1),1) = output.t+time(end);
    end
    voltage(1+size(voltage,1):size(output.voltage,1)+size(voltage,1),1) = output.voltage;
    phie(1+size(phie,1):size(output.phie,1)+size(phie,1),:) = output.phie;
    phisp(1+size(phisp,1):size(output.phisp,1)+size(phisp,1),:) = output.phisp;
    phisn(1+size(phisn,1):size(output.phisn,1)+size(phisn,1),:) = output.phisn;
    ce(1+size(ce,1):size(output.ce,1)+size(ce,1),:) = output.ce;
    csp_surf(1+size(csp_surf,1):size(output.csp_surf,1)+size(csp_surf,1),:) = output.csp_surf;
    csn_surf(1+size(csn_surf,1):size(output.csn_surf,1)+size(csn_surf,1),:) = output.csn_surf;
    SOCp(1+size(SOCp,1):size(output.soc_p_bulk,1)+size(SOCp,1),1) = output.soc_p_bulk;
    SOCn(1+size(SOCn,1):size(output.soc_n_bulk,1)+size(SOCn,1),1) = output.soc_n_bulk;
    etap(1+size(etap,1):size(output.etap,1)+size(etap,1),:) = output.etap;
    etan(1+size(etan,1):size(output.etan,1)+size(etan,1),:) = output.etan;
    Up(1+size(Up,1):size(output.Up,1)+size(Up,1),:) = output.Up;
    Un(1+size(Un,1):size(output.Un,1)+size(Un,1),:) = output.Un;
    curr_dens(1+size(curr_dens,1):size(output.curr_dens,1)+size(curr_dens,1),1) = output.curr_dens;
end

% Print performance metrics
fprintf('\nSummary of DAE solver times:\n')
fprintf('----------------------------\n')
for i = 1:cycle_num
    fprintf(['Cycle #' num2str(i) ': ' num2str(DAE_solve_time(i)) ' s\n'])
end
fprintf(['Total cycles: ' num2str(sum(DAE_solve_time)) 's\n'])
fprintf('----------------------------\n')

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
xlabel('Time [s]'); ylabel('Current, $I$ [A]');ylim([-20 20]);xlim([0 inf])
title('CC Cycling at 1C')
set(gca,'Fontsize',24);
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','latex')
set(gcf,'color','w');
set(gcf, 'Units','inches','Position', [1 1 13 8]);
legend off

% State-of-charge in positive and negative electrodes
figure; hold on; box on;
plot(time,SOCp*100,'linewidth',2)
plot(time,SOCn*100,'linewidth',2)
xlabel('Time [s]'); ylabel('State-of-charge, $SOC$ [\%]')
title('CC Cycling at 1C');legend({'$SOC_p$' '$SOC_n$'}) 
set(gca,'Fontsize',24);xlim([0 inf])
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
xlabel('Time [s]'); ylabel('[V]'); ylim([-1 1]);xlim([0 inf])
title('CC Cycling at 1C')
set(gca,'Fontsize',24);legend('location','southeast')
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','latex')
set(gcf,'color','w');
set(gcf, 'Units','inches','Position', [1 1 13 9]);

% -------------------------------------------------------------------------
% DFN_pso_0_05C.m
% -------------------------------------------------------------------------
% DFN_pso_0_05C identifies stoichiometric parameters using C/20 discharge
% data. Stoichiometric identification is commonly implemented as the first 
% step in parameter identfication to achieve accurate OCP windows for each
% electrode, which as a result determine the OCV window.
%
% NOTE: This code can be modified to identify parameters using any CC
% profile, e.g., C/5, 1C, 2C, etc. 
%
% Reference:
% [1] G. Pozzato, A. Allam, and S. Onori, “Lithium-ion battery aging dataset based on electric vehicle real-driving profiles,” Data in Brief, vol. 41, p. 107995, Apr. 2022, doi: 10.1016/j.dib.2022.107995.

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
% Option 1: Load parameters that include a previous identification 
%--------------------------------------------------------------------------
% param = load('...')
%--------------------------------------------------------------------------
% Option 2: If no previous identification, then load original nominal parameters [Chen 2020]
%--------------------------------------------------------------------------
param = Parameters_LG_INR21700_M50;

%--------------------------------------------------------------------------
% Enter mat file name where your PSO results will be stored
%--------------------------------------------------------------------------
file_name = 'pso_0_05C';

%--------------------------------------------------------------------------
% Enter names of parameters to identify (make sure names match the
% parameter names in "param" structure containing the nominal parameters)
%--------------------------------------------------------------------------
param_CC = {'theta100_p', 'theta100_n', 'theta0_p', 'theta0_n'};

%--------------------------------------------------------------------------
% Enter lower and upper bounds of parameters to identify 
%--------------------------------------------------------------------------
%   Example for parameter name "xxx":
%   lower_bounds.xxx = ...
%   upper_bounds.xxx = ...
%--------------------------------------------------------------------------
% theta100_p
lower_bounds.theta100_p = 0.22; 
upper_bounds.theta100_p = 0.34;
% theta100_n
lower_bounds.theta100_n = 0.7; 
upper_bounds.theta100_n = 1; 
% theta0_p
lower_bounds.theta0_p = 0.7; 
upper_bounds.theta0_p = 1; 
% theta0_n
lower_bounds.theta0_n = 0.015; 
upper_bounds.theta0_n = 0.04; 

%--------------------------------------------------------------------------
% Maximum particle concentration constraint
%--------------------------------------------------------------------------
% csmax_constraint = 0     -> no constraint on csp_max and csn_max 
% csmax_constraint = 1     -> constrain csp_max and csn_max 
%--------------------------------------------------------------------------
% If identifying stoichiometric parameters, you can constrain the csp_max
% and csn_max based on the theoretical capacity equation of a cell and the
% stoichiometric parameters
%--------------------------------------------------------------------------
param.csmax_constraint = 1;

%--------------------------------------------------------------------------
% PSO Settings: Enter number of particles for PSO
% (for code testing purposes, reduce particle number,e.g., particle_num = 2)
%--------------------------------------------------------------------------
particle_num = 100;

%--------------------------------------------------------------------------
% PSO Settings: Enter PSO exit condition 
% (https://www.mathworks.com/help/gads/particleswarm.html#budidgf-options)
%--------------------------------------------------------------------------
% ----------------------------------
% Option 1: Maximum iterations
% ----------------------------------
% exit_type = 'MaxIterations' -> PSO exits after maximum number of PSO iterations reached
% max_iterations = ... -> Maximum number of PSO iterations
% ----------------------------------
% Option 2: Maximum stall iterations
% ----------------------------------
% exit_type = 'MaxStallIterations' -> PSO iterations end when the relative change in best objective function value over the last MaxStallIterations iterations is less than FunctionTolerance
% max_iterations = ... -> Number of MaxStallIterations
% max_stall_tolerance = ... -> FunctionTolerance value
exit_type = 'MaxStallIterations';
max_iterations = 5;
max_stall_tolerance = 0.5e-6;

%--------------------------------------------------------------------------
% PSO Settings: Enter PSO weight 
%--------------------------------------------------------------------------
% References: 
% [1] https://www.mathworks.com/help/gads/particle-swarm-optimization-algorithm.html
% [2] G. Pozzato and S. Onori, "A General Matlab and COMSOL Co-simulation Framework for Model Parameter Optimization: Lithium-Ion Battery and Gasoline Particulate Filter Case Studies," presented at the Automotive Technical Papers, Warrendale, Pennsylvania, United States, Jul. 2023, pp. 2023-01–5047. doi: 10.4271/2023-01-5047.
%--------------------------------------------------------------------------
social_adjustment = 3.6; 
self_adjustment = 0.3;  

%--------------------------------------------------------------------------
% Load Experimental Data 
%--------------------------------------------------------------------------
%   t: Should be a vector consisting of your time experiment data      [s] (Mx1)
%   I: Should be a vector consisting of your current experiment data   [A] (Mx1) (negative current: discharging)
%   V: Should be a vector consisting of your volatge experiemntal data [V] (Mx1)
%   -> where M is the total number of data points in your experiment
%--------------------------------------------------------------------------
% C/20 capacity test conducted on LG INR21700 M50T cells
load('data_INR21700_M50T/capacity_test_data_W8_Diag1.mat')

t = t_data;
I = I_data;
V = V_data;

%--------------------------------------------------------------------------
% Choose the desired sampling time [s] to resample your experiment data 
% (the interpolated experimental data will be used in the PSO to compute the objective function)
%--------------------------------------------------------------------------
deltaT_exp = 1; 

%--------------------------------------------------------------------------
% Enter experimental data initial SOC [-] 
%--------------------------------------------------------------------------
SOC_init = 1;  

%--------------------------------------------------------------------------
% Simulation time 
%--------------------------------------------------------------------------
% NOTE: Simulation will exit automatically when the lower or upper cut-off
% voltages are reached. To simulate till the cut-off voltage condition,
% keep tf large enough, e.g., for C/20 the tf should be greater than 7200.
%--------------------------------------------------------------------------
% Initial simulation time [s]
t0 = 0;
% Final simulation time [s]
tf = 1e6;

%--------------------------------------------------------------------------
% Print cost function value (useful if you want to see the J value calculated for every particle)
%--------------------------------------------------------------------------
% param.J_print = 0     -> supress J value during PSO
% param.J_print = 1     -> print J value during PSO
%--------------------------------------------------------------------------
param.J_print = 1;

%--------------------------------------------------------------------------
% Current mode (do not change)
%--------------------------------------------------------------------------
% param.CurrentMode = 1     -> constant current (CC)
% param.CurrentMode = 2     -> variable current, e.g., UDDS
%--------------------------------------------------------------------------
param.CurrentMode = 1;

%% Processing experimental data for PSO

% Interpolate experimental data according to user defined sampling time
t_exp = (0:deltaT_exp:t(end))'; 
I_exp = interp1(t,I,t_exp);     % Negative current: discharging
V_exp = interp1(t,V,t_exp);

% Calculate discharge capacity during experimental data        
Q_dis_exp = abs(trapz(t_exp/3600,-I_exp));
    
% Calculate SOC using coulomb counting
SOC_exp = SOC_init - cumtrapz(t_exp/3600,-I_exp)/param.Q_nom;

% Ensure that SOC is not below 0% or above 100% 
SOC_exp(SOC_exp<0) = 0;
SOC_exp(SOC_exp>1) = 1;

% Define SOCp and SOCn
SOCp_exp = SOC_exp;
SOCn_exp = SOC_exp;

%% Processing simulation settings

current = mean(I_exp);
param.currentDensity = current/param.Acell;
param.SOC_init = SOC_init;
param.t0 = t0;
param.tf = tf;

%% Processing parameters for PSO

init_pos = zeros(1,length(param_CC));
lb = zeros(1,length(param_CC));
ub = zeros(1,length(param_CC));
% Set initial position of parameters as the nominal values
for i = 1:length(param_CC)
    init_pos(i) = param.(param_CC{i});
    lb(i) = lower_bounds.(param_CC{i});
    ub(i) = upper_bounds.(param_CC{i});
end

%% Print initial position and bounds
fprintf('Displaying parameter settings...\n')
for i = 1:length(init_pos)
    fprintf('------------------------\n')
    fprintf([param_CC{i} ':\n'])
    fprintf([num2str(lb(i)) '(lower) | ' num2str(init_pos(i)) '(initial) | ' num2str(ub(i)) '(upper)\n'])
end

%% PSO setting 
% Make sure to have "UseVectorized = true" to run PSO simulations in parallel
% for faster computation time (requires Parallel Computing Toolbox)

% FunctionTolerance: Iterations end when the relative change in best 
% objective function value over the last MaxStallIterations iterations is 
% less than options.FunctionTolerance.max
options = optimoptions('particleswarm',...
                       'UseVectorized',true,...
                       'SwarmSize',particle_num,...
                       'Display','iter',...
                       exit_type,max_iterations,...
                       'FunctionTolerance',max_stall_tolerance,...
                       'InitialSwarmMatrix',init_pos,...
                       'SocialAdjustmentWeight',social_adjustment,... 
                       'SelfAdjustmentWeight',self_adjustment); 

%% RUN PSO 
fprintf('\n\nStart running PSO!...\n')

% To pass additional parameters to PSO objective function, use anonymous function method
fcn = @(x)DFN_obj_CC(x,param,param_CC,t_exp,I_exp,V_exp,SOCp_exp,SOCn_exp,Q_dis_exp);

% fval: objective function value
% exitflag: exit condition stopping condition
% output: information about optimization
[x_opt_CC,fval,exitflag,output] = particleswarm(fcn,length(init_pos),lb,ub,options);

save(['DFN_identification_results/' file_name '.mat'],'x_opt_CC','lb','ub','init_pos','param_CC','Q_dis_exp','options','fval','exitflag','output')

%% Print identified values
fprintf('Displaying identified values...\n')
for i = 1:length(param_CC)
    fprintf('------------------------\n')
    fprintf([param_CC{i} ':\n'])
    fprintf(['Identified value: ' num2str(x_opt_CC(i)) '\n'])
    fprintf([num2str(lb(i)) '(lower) | ' num2str(init_pos(i)) '(initial) | ' num2str(ub(i)) '(upper)\n'])
end

%% Simulate using identified parameters

% Update parameters to identified parameters
for i = 1:length(param_CC)
    param.(param_CC{i}) = x_opt_CC(i);
end
if param.csmax_constraint == 1
    fprintf('\nStoichiometric identification: maximum solid concentration constraint was enabled\n')
    param.ctp = 3600*param.Q_nom/(param.ep_s*param.F*param.lp*param.Acell*(param.theta0_p-param.theta100_p));
    param.ctn = 3600*param.Q_nom/(param.en_s*param.F*param.ln*param.Acell*(param.theta100_n-param.theta0_n));
    fprintf(['csp_max = ' num2str(param.ctp) '\n'])
    fprintf(['csn_max = ' num2str(param.ctn) '\n'])
end

% Do not change code here (updating constrained parameters)
%-----------------------------------------------------------------------
param.ep_s=1-param.ep;
param.en_s=1-param.en;
param.ap=(3/param.Rpp)*param.ep_s;
param.an=(3/param.Rpn)*param.en_s;
param.currentDensity = current/param.Acell;
%-----------------------------------------------------------------------

% Save updated parameters based on identified parameters in mat file
save(['DFN_identification_results/' file_name '.mat'],'param','t_exp','V_exp','I_exp','SOCp_exp','SOCn_exp','SOC_init','t0','tf','-append')

x_init = [];
results = runModel(x_init,param,t0,tf,SOC_init,param.currentDensity);

%% Visualization: plot identified results against experimental data

% Simulated data
t_sim = results.t;
V_sim = results.voltage;
SOCp_sim = results.soc_p_bulk;
SOCn_sim = results.soc_n_bulk;

figure; hold on; box on
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
fprintf('Displaying objective function values...\n')
fprintf('------------------------\n')
fprintf(['J_V =' num2str(J1) ' [-]\n'])
fprintf(['J_V =' num2str(J1_mV*1000) ' [mV]\n'])
fprintf(['J_SOCp =' num2str(J2*100) ' [%%]\n'])
fprintf(['J_SOCn =' num2str(J3*100) ' [%%]\n'])
fprintf(['J_tot =' num2str(J1+J2+J3) ' [-]\n'])

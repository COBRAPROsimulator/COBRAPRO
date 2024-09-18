% -------------------------------------------------------------------------
% DFN_pso_HPPC.m
% -------------------------------------------------------------------------
% DFN_pso_HPPC identifies parameters using experimentally obtained 
% HPPC voltage profile. The C/20 identification results from "DFN_pos_0_05C.m"
% are utilized. All other parameters are obtained from [Chen 2020].
%
% Reference:
% [1] G. Pozzato, A. Allam, and S. Onori, “Lithium-ion battery aging dataset based on electric vehicle real-driving profiles,” Data in Brief, vol. 41, p. 107995, Apr. 2022, doi: 10.1016/j.dib.2022.107995.
% [2] C.-H. Chen, F. Brosa Planella, K. O’Regan, D. Gastol, W. D. Widanage, and E. Kendrick, “Development of Experimental Techniques for Parameterization of Multi-scale Lithium-ion Battery Models,” J. Electrochem. Soc., vol. 167, no. 8, p. 080534, Jan. 2020, doi: 10.1149/1945-7111/ab9050.

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
load('identified_parameters_0_05C.mat','param')
%--------------------------------------------------------------------------
% Option 2: If no previous identification, then load original nominal parameters [Chen 2020]
%--------------------------------------------------------------------------
% param = Parameters_LG_INR21700_M50;

%--------------------------------------------------------------------------
% Enter mat file name where your PSO results will be stored
%--------------------------------------------------------------------------
file_name = 'pso_HPPC';

%--------------------------------------------------------------------------
% Option 1: Enter names of parameters to identify (make sure names match the
% parameter names in "param" structure containing nominal parameters)
%--------------------------------------------------------------------------
% param_HPPC = {'Dsp' 'Dsn' 'kp' 'kn' 'De' 'Kappa'};
%--------------------------------------------------------------------------
% Option 2: Load identifiable parameters from identifiability analysis
% conducted in "Examples/Parameter_Identifiability_Analysis/DFN_LSA_Corr_HPPC.m"
%--------------------------------------------------------------------------
load('HPPC_identifiable_params.mat')
% Enter desired beta value
beta_value = 0.95;
% Load results from LSA and correlation analysis 
if ismember(beta_value,beta_corr)
    corr_ind = beta_value == beta_corr;
else
    fprintf('Invalid correlation coefficient threshold.\nPlease choose a correlation threshold that was tested in your parameter identifiability analysis.\n')
    return
end
param_HPPC = corr_identifiable_vec{corr_ind};

%--------------------------------------------------------------------------
% Enter lower and upper bounds of parameters to identify 
%--------------------------------------------------------------------------
%   Example for parameter name "xxx":
%   lower_bounds.xxx = ...
%   upper_bounds.xxx = ...
%--------------------------------------------------------------------------
% Dsp
pct = 0.2; % perturbation coeff
lower_bounds.Dsp = 10^(log10(param.Dsp)*(1+pct));
upper_bounds.Dsp = 10^(log10(param.Dsp)*(1-pct));
% kn
pct = 0.3; % perturbation coeff
lower_bounds.kn = 10^(log10(param.kn)*(1+pct));
upper_bounds.kn = 10^(log10(param.kn)*(1-pct));
% c0
lower_bounds.c0 = 500;
upper_bounds.c0 = 1500;
% Dsn
pct = 0.2; % perturbation coeff
lower_bounds.Dsn = 10^(log10(param.Dsn)*(1+pct));
upper_bounds.Dsn = 10^(log10(param.Dsn)*(1-pct));

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
% HPPC test conducted on LG INR21700 M50T cells
load('data_INR21700_M50T/HPPC_data_W8_Diag1.mat')

t = t_data;
I = I_data;
V = V_data;

%--------------------------------------------------------------------------
% Enter initial SOC of your HPPC data
%--------------------------------------------------------------------------
SOC_init = 1;   % [-]      

%--------------------------------------------------------------------------
% Choose the desired sampling time [s] to resample your experiment data 
% (the interpolated experimental data will be used in the PSO to compute the objective function)
%--------------------------------------------------------------------------
deltaT_exp = 1; % [s]

%--------------------------------------------------------------------------
% Load HPPC profile input current vector [A] and simulation time [s] 
%--------------------------------------------------------------------------
% HPPC consists of simulating current pulses and rest periods during CC discharging.
% Here we simulate HPPC by concatenating multiple CC simulations, called CC segments.
%--------------------------------------------------------------------------
% time_vec: Should be a vector where each entry defines the time period for each CC segment [s] (Px1)
% curr_vec: Should be a vector where each entry defines the current for that CC segment [A] (Px1) (negative current: discharging)
% -> where P is the number of CC segments for your HPPC experiment
% For example, since HPPC starts with a 1C discharge for 360 s followed by
% a rest period for 3600 s:
%   -> time_vec(1) = 360
%   -> curr_vec(1) = -4.8504
%   -> time_vec(2) = 3600
%   -> curr_vec(2) = 0
%--------------------------------------------------------------------------
% Load HPPC_W8_Diag1.mat file that consists of curr_vec and time_vec
load('HPPC_W8_Diag1.mat') % Data extrated from [1] 

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
% param.CurrentMode = 1     -> constant current
% param.CurrentMode = 2     -> variable current, e.g., UDDS
%--------------------------------------------------------------------------
param.CurrentMode = 1;

%% Processing experimental data for PSO

% Interpolate experimental data according to user defined sampling time
t_exp = (0:deltaT_exp:t(end))'; 
I_exp = interp1(t,I,t_exp);     % Negative current: discharging
V_exp = interp1(t,V,t_exp);

% Calculate discharge capacity during HPPC experimental data        
Q_dis_exp_HPPC = abs(trapz(t_exp/3600,-I_exp));

% Calculate SOC using coulomb counting (using nominal capacity or Q_dis)
% For more information, refer to:
% [2] S. Ha, G. Pozzato, and S. Onori, "Electrochemical characterization tools for lithium-ion batteries," J Solid State Electrochem, Nov. 2023, doi: 10.1007/s10008-023-05717-1.
SOC_exp = SOC_init - cumtrapz(t_exp/3600,-I_exp)/param.Q_nom;

% Ensure that SOC is not below 0% or above 100% 
SOC_exp(SOC_exp<0) = 0;
SOC_exp(SOC_exp>1) = 1;

% Define SOCp and SOCn
SOCp_exp = SOC_exp;
SOCn_exp = SOC_exp;

%% HPPC simulation

% Define current density vector for each CC segment
currentDensity_vec = curr_vec/param.Acell;
param.curr_dens_vec = curr_vec/param.Acell;

% Number of CC segments (equal to P)
curr_segments = length(curr_vec); 

% Initial simulation time [s]
t0 = 0;
% Final simulation time vector [s]
tf_vec = time_vec; 

%% Processing HPPC parameters for PSO

init_pos = zeros(1,length(param_HPPC));
lb = zeros(1,length(param_HPPC));
ub = zeros(1,length(param_HPPC));
% Set initial position of parameters as the nominal values
for i = 1:length(param_HPPC)
    init_pos(i) = param.(param_HPPC{i});
    lb(i) = lower_bounds.(param_HPPC{i});
    ub(i) = upper_bounds.(param_HPPC{i});
end

%% Print initial position and bounds
fprintf('\nDisplaying parameter settings...\n')
for i = 1:length(init_pos)
    fprintf('------------------------\n')
    fprintf([param_HPPC{i} ':\n'])
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
fprintf('\nStart running PSO!...\n')

% To pass additional parameters to PSO objective function, use anonymous function method
fcn = @(x)DFN_obj_HPPC(x,param,param_HPPC,SOC_init,t0,tf_vec,curr_segments,t_exp,curr_vec,V_exp,SOCp_exp,SOCn_exp,Q_dis_exp_HPPC);

% fval: objective function value
% exitflag: exit condition stopping condition
% output: information about optimization
[x_opt_HPPC,fval,exitflag,output] = particleswarm(fcn,length(init_pos),lb,ub,options);

save(['Examples/Parameter_Identification_Results/' file_name '.mat'],'x_opt_HPPC','lb','ub','init_pos','param_HPPC','param','Q_dis_exp_HPPC')

%% Print identified values
fprintf('\nDisplaying identified values...\n')
for i = 1:length(param_HPPC)
    fprintf('------------------------\n')
    fprintf([param_HPPC{i} ':\n'])
    fprintf(['Identified value: ' num2str(x_opt_HPPC(i)) '\n'])
    fprintf([num2str(lb(i)) '(lower) | ' num2str(init_pos(i)) '(initial) | ' num2str(ub(i)) '(upper)\n'])
end

%% Simulate using identified parameters

% Update parameters to identified parameters
for i = 1:length(param_HPPC)
    param.(param_HPPC{i}) = x_opt_HPPC(i);
end

% Do not change code here (updating constrained parameters)
%-----------------------------------------------------------------------
param.ep_s=1-param.ep;
param.en_s=1-param.en;
param.ap=(3/param.Rpp)*param.ep_s;
param.an=(3/param.Rpn)*param.en_s;
param.curr_dens_vec = curr_vec/param.Acell;
%-----------------------------------------------------------------------

% Save updated parameters based on identified parameters in mat file
save(['Examples/Parameter_Identification_Results/' file_name '.mat'],'param','t_exp','V_exp','I_exp','SOCp_exp','SOCn_exp','SOC_init','t0','tf_vec','currentDensity_vec','curr_segments','-append')

fprintf('\nSimulating model using identified parameters...\n')
[t_sim,V_sim,SOCp_sim,SOCn_sim,curr_dens_sim,...
    delta_phie_sim,delta_eta_sim,delta_OCV_sim]=HPPC_sim(param,SOC_init,t0,tf_vec,currentDensity_vec,curr_segments);

%% Visualization: plot identified results against experimental data

myexp_color = [0.5138 0.7332 0.8587];

figure; hold on; box on; set(gcf,'color','w');
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

%% Uncomment to plot internal states of simulated HPPC results (helpful to better understand overpotentials contributing to voltage)
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

% figure; hold on
% % Plot simulated voltage and simulated OCV 
% yyaxis left
% plot(t_sim/3600,V_sim,'Linewidth',2,'DisplayName','$V_{cell}^{sim}=\Delta OCP+\Delta \phi_e+\Delta \eta + IR_c$')
% plot(t_sim/3600,delta_OCV_sim,'Linewidth',2,'DisplayName','$\Delta OCV$')
% ylabel('Voltage [V]'); xlabel('Time [h]')
% % Plot delta_eta, delta_phie, delta_OCV, annd Rc*I
% yyaxis right
% plot(t_sim/3600,delta_phie_sim,'Linewidth',2,'DisplayName','$\Delta\phi_{e}$')
% plot(t_sim/3600,delta_eta_sim,'Linewidth',2,'DisplayName','$\Delta\eta$')
% plot(t_sim/3600,curr_dens_sim*param.Rc*param.Acell,'Linewidth',2,'DisplayName','$IR_c$')
% ylabel('[V]');
% xlim([0 inf]);
% yyaxis left
% ylim([2.5 4.2])
% yyaxis right
% ylim([-0.2 0.2])
% title('Simulated Voltage and Internal Variables')
% legend('location','southwest')
% set(gcf,'color','w');set(gca,'Fontsize',32);
% set(findall(gcf,'-property','interpreter'),'interpreter','latex')
% set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','latex')
% set(gcf, 'Units','inches','Position', [1 1 13 8]);

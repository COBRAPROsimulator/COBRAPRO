% -------------------------------------------------------------------------
% test_2_comsolValidation.m
% -------------------------------------------------------------------------
%   test_2_comsolValidation checks COBRAPRO's output against COMSOL for 1C
%   discharge. This script automatically checks the RMSE voltage against COMSOL and 
%   tells the user if the test was successful. A figure comparing 
%   COBRAPRO and COMSOL is provided as an additional reference.

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

close all;clc;clear;

% Load COMSOL 1C discharge data (same parameters used as COBRAPRO results)
load('COMSOL_1C.mat')

% Initialize all model parameters
param = Parameters_LG_INR21700_M50_comsolValidation;

% Initial conditions
x_init = [];

% Define current density [A.m-2]
Crate = -1; 
current = Crate*param.Q_nom; 
currentDensity = current/param.Acell;

% Initial SOC [-]
SOC_init = 1;

% Simulation time [s]
t0 = 0; % Initial simulation time [s]
tf = 4000; % Final simulation time [s]

% Run COBRAPRO
[output,param] = runModel(x_init,param,t0,tf,SOC_init,currentDensity);

figure;hold on;
plot(comsol_time,comsol_voltage,'LineWidth',2,'DisplayName','COMSOL')
plot(output.t,output.voltage,'--','LineWidth',2,'DisplayName','COBRAPRO')
xlabel('time [s]');ylabel('Voltage [V]');legend('Location','northeast')
title('COMSOL Validation at 1C Discharge')
set(gca,'Fontsize',32);set(gcf,'color','w');
set(findall(gcf,'-property','interpreter'),'interpreter','latex')
set(findall(gcf,'-property','ticklabelinterpreter'),'ticklabelinterpreter','latex')
set(gcf, 'Units','inches','Position', [1 1 13 8]);

% RMSE Error calculation
cobrapro_voltage = interp1(output.t,output.voltage,output.t);
min_len = min(length(cobrapro_voltage),length(comsol_voltage));
error = rms(cobrapro_voltage(1:min_len)-comsol_voltage(1:min_len));

error_test = 1e-2;
if error < error_test
    fprintf('\ntest_2 successful: COBRAPRO is working as expected! Results validated against COMSOL.\n\n');
else 
    fprintf('\ntest_2 unsuccessful: COBRAPRO results not validated against COMSOL.\n\n');
end

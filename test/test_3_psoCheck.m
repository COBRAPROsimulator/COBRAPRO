% -------------------------------------------------------------------------
% test_3_psoCheck.m
% -------------------------------------------------------------------------
%   test_3_psoCheck checks that the MATLAB Global Optimization Toolbox and
%   MATLAB Parallel Computing Toolbox are installed, which are required to
%   run the parameter optimization routine scripts in the folder
%   Examples/Parameter_Identification_Routines. This test also checks that
%   the particle swarm optimization (PSO) in parallel works as expected. 

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

% Check that MATLAB Global Optimization Toolbox is installed
GADS_check = license('test','GADS_Toolbox');
if GADS_check == 1
    fprintf('Global Optimization Toolbox is installed!\n\n')
else
    fprintf(['Global Optimization Toolbox is not installed...\n' ...
        'Please install Global Optimization Toolbox to use COBRAPRO and complete test_3_psoCheck.\n']);
    return
end

% Check that MATLAB Parallel Computing Toolbox is installed
GADS_check = license('test','Distrib_Computing_Toolbox');
if GADS_check == 1
    fprintf('Parallel Computing Toolbox is installed!\n\n')
else
    fprintf(['Parallel Computing Toolbox is not installed...\n' ...
        'Please install Parallel Computing Toolbox to use COBRAPRO and complete test_3_psoCheck.\n']);
    return
end

% PSO settings
init_pos = [0.5 0.6];
lb = [0.1 0.2];
ub = [0.8 0.9];
particle_num = 3;
options = optimoptions('particleswarm',...
                       'UseVectorized',true,...
                       'SwarmSize',particle_num,...
                       'Display','iter',...
                       'MaxIterations',5,...
                       'InitialSwarmMatrix',init_pos);  

% RUN PSO 
fprintf('Start running PSO...\n')

% To pass additional parameters to PSO objective function, use anonymous function method
number = 2;
fcn = @(x)PSO_obj_test(x,number);

% fval: objective function value
% exitflag: exit condition stopping condition
% output: information about optimization
[x_opt,fval,exitflag,output] = particleswarm(fcn,length(init_pos),lb,ub,options);

if isnumeric(x_opt) && (length(x_opt) == length(init_pos))
    fprintf('\ntest_3 successful: PSO in parallel is working as expected!\n\n');
else 
    fprintf('\ntest_3 unsuccessful: PSO in parallel is not working.\n\n');
end

function J = PSO_obj_test(x,pso_workers)
%   DFN_obj_CC is the objective function that outputs the RMSE between
%   experimental and simulated voltage and SOC for constant current (CC)
%   experiments
%
%   Inputs:
%       x: PSO particle value matrix
%           -> size nxm, where n is the number of particles, m is the number of parameters to identify
%       number: number of 
%
%   Output:
%       J: objective function value, calculated as J = 0.5 - (a+b)

% Preallocate J_tot for memory efficiency
J = zeros(size(x,1),1);

% Extract parameters from particles
a = x(:,1);
b = x(:,2);

parfor (i = 1:size(x,1), pso_workers)
    J(i) = 0.5 - (a(i)+b(i));
end

end
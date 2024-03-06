function J_tot = DFN_obj_CC(x,param,param_0_05C,t_exp,I_exp,V_exp,SOCp_exp,SOCn_exp,Q_dis_exp)
%   DFN_obj_CC is the objective function that outputs the RMSE between
%   experimental and simulated voltage and SOC for constant current (CC)
%   experiments
%
%   Inputs:
%       x: PSO particle value matrix
%           -> size nxm, where n is the number of particles, m is the number of parameters to identify
%       param: nominal parameter structure
%       param_0_05C: array containing names of parameters to be identified
%       t_exp: experimental time vector    (Mx1)
%       I_exp: experimental current vector (Mx1)
%       V_exp: experimental current vector (Mx1)
%       SOCp_exp: experimental SOCp vector (Mx1)
%       SOCn_exp: experimental SOCn vector (Mx1)
%           -> where M is the total number of data points in your experiment
%       Q_dis_exp: Experimental discharge capacity [Ah]
%
%   Output:
%       J_tot: objective function value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COBRAPRO: Co-simulation Battery Modeling for Accelerated Parameter Optimization

% Copyright (c) 2024 CO-simulation BatteRy modeling for Accelerated PaRameter Optimization (COBRAPRO)
% COBRAPRO is freely distributed under the MIT License 
% v1.0.0.: Released March 1, 2024 

% Main contributors: 
% Sara Ha (sungyeon.sara.ha@stanford.edu)
% Simona Onori (sonori@stanford.edu)
% Stanford Energy Control Group, Stanford University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Experimental/simulation current (only valued for CC)
current = mean(I_exp);

% Preallocate param_copies for memory efficiency
param_copies(size(x,1)) = param;

% Create copies of param with updated parameter values from the particles
parfor (i = 1:size(x,1), param.pso_workers)
% for i = 1:size(x,1)
    param_copies(i) = param;
    % Delegate particles to parameters
    for j = 1:length(param_0_05C)
        param_copies(i).(param_0_05C{j}) = x(i,j);
    end
    if param_copies(i).csmax_constraint == 1
        param_copies(i).ctp = 3600*param_copies(i).Q_nom/(param_copies(i).ep_s*param_copies(i).F*param_copies(i).lp*param_copies(i).Acell*(param_copies(i).theta0_p-param_copies(i).theta100_p));
        param_copies(i).ctn = 3600*param_copies(i).Q_nom/(param_copies(i).en_s*param_copies(i).F*param_copies(i).ln*param_copies(i).Acell*(param_copies(i).theta100_n-param_copies(i).theta0_n));
    end
    % Do not change code here (updating constrained parameters)
    %----------------------------------------------------------------------
    param_copies(i).ep_s=1-param_copies(i).ep;
    param_copies(i).en_s=1-param_copies(i).en;
    param_copies(i).ap=(3/param_copies(i).Rpp)*param_copies(i).ep_s;
    param_copies(i).an=(3/param_copies(i).Rpn)*param_copies(i).en_s;
    param_copies(i).currentDensity = current/param_copies(i).Acell;
    %----------------------------------------------------------------------
end

% Preallocate J_tot for memory efficiency
J_tot = zeros(size(x,1),1);

x_init=[];

parfor (i = 1:size(x,1), param.pso_workers)
% for i = 1:size(x,1)
    temp_params = param_copies(i);
    try
        output = runModel(x_init,temp_params,temp_params.t0,temp_params.tf,temp_params.SOC_init,temp_params.currentDensity);
        % This means that the model failed to initialize algebraic variables
        if isnan(output.t)
            J_tot(i) = 40;
        else
            % Simulated output
            t_sim = output.t;
            V_sim = output.voltage;
            SOCp_sim = output.soc_p_bulk;
            SOCn_sim = output.soc_n_bulk;
            % Calculate discharge capacity of simulated data
            Q_dis_sim = trapz(t_sim/3600,abs(current)*ones(length(t_sim),1));
            % Check constraint on discharge capacity
            if abs((Q_dis_sim-Q_dis_exp))/Q_dis_exp > 0.01
                J_tot(i) = 10;
            else
                % Interpolate sim output to match the timestep of voltage experiment
                V_sim_obj = interp1(t_sim, V_sim, t_exp);
                SOCp_sim_obj = interp1(t_sim, SOCp_sim, t_exp);
                SOCn_sim_obj = interp1(t_sim, SOCn_sim, t_exp);
                % Remove NaN from interpolated data
                V_sim_obj=V_sim_obj(~isnan(V_sim_obj));
                SOCp_sim_obj=SOCp_sim_obj(~isnan(SOCp_sim_obj));
                SOCn_sim_obj=SOCn_sim_obj(~isnan(SOCn_sim_obj));
                % Use shorter length to compare
                len1=length(V_sim_obj);len2=length(V_exp);
                len_V=min(len1,len2);
                len1=length(SOCp_sim_obj);len2=length(SOCp_exp);
                len_SOCp=min(len1,len2);
                len1=length(SOCn_sim_obj);len2=length(SOCn_exp);
                len_SOCn=min(len1,len2);
                % Objective function
                J1 = rms((V_exp(1:len_V)-V_sim_obj(1:len_V))./V_exp(1:len_V));
                J2 = rms(SOCp_exp(1:len_SOCp)-SOCp_sim_obj(1:len_SOCp));
                J3 = rms(SOCn_exp(1:len_SOCn)-SOCn_sim_obj(1:len_SOCn));
                J_tot(i)=J1+J2+J3;
            end
        end
    catch 
        J_tot(i)=50;  
    end
    if temp_params.J_print == 1
        fprintf(['J_tot ',num2str(J_tot(i)),'\n']);   
    end
end
end

function J_tot = DFN_obj_HPPC(x,param,param_HPPC,SOC_init,t0,tf_vec,curr_segments,t_exp,curr_vec,V_exp,SOCp_exp,SOCn_exp,Q_dis_exp_HPPC)
%   DFN_obj_HPPC is the objective function that outputs the RMSE between
%   experimental and simulated voltage and SOC
%
%   Inputs:
%       x: PSO particle value matrix
%           -> size nxm, where n is the number of particles, m is the number of parameters to identify
%       param: nominal parameter structure
%       param_HPPC: array containing names of parameters to be identified
%       SOC_init: initial SOC [-]
%       t0: initial time [s]
%       tf_vec: final time vector for each CC segment [s]
%       curr_segments: number of CC segments (let's say is equal to Q)
%       t_exp: experimental time vector    [s] (Mx1)
%       curr_vec: vector containing current value for each CC segment in HPPC [A] (1xQ)
%       V_exp: experimental current vector [V] (Mx1)
%       SOCp_exp: experimental SOCp vector [-] (Mx1)
%       SOCn_exp: experimental SOCn vector [-] (Mx1)
%           -> where M is the total number of data points in your experiment
%       Q_dis_exp_HPPC: Experimental discharge capacity [Ah]
%
%   Output:
%       J_tot: objective function value [-]

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

% Preallocate param_copies for memory efficiency
param_copies(size(x,1)) = param;

% Create copies of param with updated parameter values from the particles
parfor (i = 1:size(x,1), param.pso_workers)
% for i = 1:size(x,1)
    param_copies(i) = param;
    % Delegate particles to parameters
    for j = 1:length(param_HPPC)
        param_copies(i).(param_HPPC{j}) = x(i,j);
    end
    % Do not change code here (updating constrained parameters)
    %----------------------------------------------------------------------
    param_copies(i).ep_s=1-param_copies(i).ep;
    param_copies(i).en_s=1-param_copies(i).en;
    param_copies(i).ap=(3/param_copies(i).Rpp)*param_copies(i).ep_s;
    param_copies(i).an=(3/param_copies(i).Rpn)*param_copies(i).en_s;
    param_copies(i).curr_dens_vec = curr_vec./param_copies(i).Acell;
    %----------------------------------------------------------------------
end

% Preallocate J_tot for memory efficiency
J_tot = zeros(size(x,1),1);

parfor (i = 1:size(x,1),param.pso_workers)
% for i = 1:size(x,1)
    temp_params = param_copies(i);
    try
        %----------------------------------------------------------------------
        % 1. Start HPPC simulation
        %----------------------------------------------------------------------
        [t_sim,V_sim,SOCp_sim,SOCn_sim,curr_dens_sim,~,~,~,init_fail]=HPPC_sim(temp_params,SOC_init,t0,tf_vec,temp_params.curr_dens_vec,curr_segments);

        %----------------------------------------------------------------------
        % 2. Define objective function
        %----------------------------------------------------------------------
        %----------------------------------------------------------------------
        % If the HPPC simulation could not complete due to initialization issues during one of the CC segment simulations
        %----------------------------------------------------------------------
        if init_fail ~= 0
            J_tot(i) = 40;
            message = 'Simulation could not complete due to initialization issues during HPPC';
        %----------------------------------------------------------------------
        % If HPPC simulation did not exit due to CC segment initialization issues
        %----------------------------------------------------------------------
        else
            % Calculate discharge capacity during HPPC simulation
            Q_dis_sim_HPPC = abs(trapz(t_sim/3600,-curr_dens_sim*temp_params.Acell));
            %----------------------------------------------------------------------
            % If simulation ended preemptively due to exit conditions 
            %----------------------------------------------------------------------
            if abs((Q_dis_sim_HPPC - Q_dis_exp_HPPC))/Q_dis_exp_HPPC > 0.01
                % If simulation ended preemptively due to exit conditions being
                % triggered, then the discharge capacity of the simulated HPPC
                % will be less than that of the experimental HPPC
                J_tot(i) = 30;
                message = 'Simulation could not complete since simulation and experiment discharge capacities do not match.';
            %----------------------------------------------------------------------
            % If there was an error in calculating the Q_dis_sim_HPPC
            %----------------------------------------------------------------------
            elseif isnan(Q_dis_sim_HPPC)
                J_tot(i) = 30;
                message = 'Simulation capacity is NaN.';
            %----------------------------------------------------------------------
            % If all CC segments successfully simulated without triggering any exit conditions
            %----------------------------------------------------------------------
            else
                % Interpolate sim output to match the timestep of voltage experiment
                V_sim_obj = interp1(t_sim, V_sim, t_exp);
                SOCp_sim_obj = interp1(t_sim, SOCp_sim, t_exp);
                SOCn_sim_obj = interp1(t_sim, SOCn_sim, t_exp);
                % Remove NaN from interpolated data
                V_sim_obj=V_sim_obj(~isnan(V_sim_obj));
                SOCp_sim_obj=SOCp_sim_obj(~isnan(SOCp_sim_obj));
                SOCn_sim_obj=SOCn_sim_obj(~isnan(SOCn_sim_obj));
                % use shorter length to compare
                len1=length(V_sim_obj);
                len2=length(V_exp);
                len=min(len1,len2);
                len1=length(SOCp_sim_obj);len2=length(SOCp_exp);
                len_SOCp=min(len1,len2);
                len1=length(SOCn_sim_obj);len2=length(SOCn_exp);
                len_SOCn=min(len1,len2);
                % Objective function
                J1 = rms((V_exp(1:len)-V_sim_obj(1:len))./V_exp(1:len));
                J2 = rms(SOCp_sim_obj(1:len_SOCp)-SOCp_exp(1:len_SOCp));
                J3 = rms(SOCn_sim_obj(1:len_SOCn)-SOCn_exp(1:len_SOCn));
                J_tot(i) = J1+J2+J3;
                message = 'No error.';
            end
        end
    catch 
        J_tot(i)=50;  
        message = 'Parameters resulted in simulation error.';
    end
    if temp_params.J_print == 1
        fprintf(['J_tot ',num2str(J_tot(i)), '. Error: ' message '\n']);     
    end
end

end
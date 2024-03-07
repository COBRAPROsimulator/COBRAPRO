function [t_nom,V_nom,SOCp_nom,SOCn_nom,...
    t_l_all,V_l_all,SOCp_l_all,SOCn_l_all,...
    t_u_all,V_u_all,SOCp_u_all,SOCn_u_all,...
    S_norm_V, S_norm_SOCp, S_norm_SOCn, S_norm_V_SOC]...
    = DFN_LSA_sim_HPPC(param,param_LSA_HPPC,theta_nom,pct)
%   DFN_LSA_sim_HPPC conducts local sensitivity analysis (LSA) by locally 
%   perturbing the parameters of interest around a nominal value for HPPC simulations. 
%
%   Inputs:
%       param: parameter structure
%       param_LSA_HPPC: cell consisting of parameter names to conduct LSA on
%       theta_nom: vector consisting of nominal values of the parameters to conduct LSA on
%       pct: perturbation coefficient [-]
%
%   Output:
%       t_nom: time vector simulated using nominal parameters [s]
%       V_nom: voltage vector simulated using nominal parameters [V]
%       SOCp_nom: Positive electrode SOC vector simulated using nominal parameters [-]
%       SOCn_nom: Negative electrode SOC vector simulated using nominal parameters [-]
%       t_l_all: cell array where each cell contains the time vector simulated using the lower perturbed parameter [s]
%       V_l_all: cell array where each cell contains the voltage vector simulated using the lower perturbed parameter [V]
%       SOCp_l_all: cell array where each cell contains the positive electrode SOC vector simulated using the lower perturbed parameter [-]
%       SOCn_l_all: cell array where each cell contains the negative electrode SOC vector simulated using the lower perturbed parameter [-]
%       t_u_all: cell array where each cell contains the time vector simulated using the upper perturbed parameter [s]
%       V_u_all: cell array where each cell contains the voltage vector simulated using the upper perturbed parameter [V]
%       SOCp_u_all: cell array where each cell contains the positive electrode SOC vector simulated using the upper perturbed parameter [-]
%       SOCn_u_all: cell array where each cell contains the negative electrode SOC vector simulated using the upper perturbed parameter [-]
%       S_norm_V: normalized voltage sensitivity vector [-]
%       S_norm_SOCp: normalized positive electrode SOC sensitivity vector [-]
%       S_norm_SOCn: normalized negative electrode SOC sensitivity vector [-]
%       S_norm_V_SOC: normalized total sensitivity vector (including voltage, positive electrode SOC, and negative electrode SOC sensitivities) [-]

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

%----------------------------------------------------------------------
% Simulate using nominal parameters
%----------------------------------------------------------------------
fprintf('Simulating with nominal parameters...\n')
[t_nom,V_nom,SOCp_nom,SOCn_nom,~,~,~,~]=HPPC_sim(param,param.SOC_init,param.t0,param.tf_vec,param.curr_dens_vec,param.curr_segments);
fprintf('Done.\n')

%----------------------------------------------------------------------
% Simulate HPPC LSA
%----------------------------------------------------------------------
% Preallocate cells for memory efficiency
S_norm_V=cell(1,length(param_LSA_HPPC));
S_norm_SOCp=cell(1,length(param_LSA_HPPC));
S_norm_SOCn=cell(1,length(param_LSA_HPPC));
S_norm_V_SOC=cell(1,length(param_LSA_HPPC));
t_l_all=cell(1,length(param_LSA_HPPC));
V_l_all=cell(1,length(param_LSA_HPPC));
SOCp_l_all=cell(1,length(param_LSA_HPPC));
SOCn_l_all=cell(1,length(param_LSA_HPPC));
t_u_all=cell(1,length(param_LSA_HPPC));
V_u_all=cell(1,length(param_LSA_HPPC));
SOCp_u_all=cell(1,length(param_LSA_HPPC));
SOCn_u_all=cell(1,length(param_LSA_HPPC));

% Preallocate param_copies for memory efficiency
param_copies_l(length(param_LSA_HPPC)) = param;
param_copies_u(length(param_LSA_HPPC)) = param;

% Create copies of param with updated perturbed parameter values
parfor (i = 1:length(param_LSA_HPPC), param.pso_workers)
    %----------------------------------------------------------------------
    % Update param structures with the lower bound perturbed parameters
    %----------------------------------------------------------------------
    param_copies_l(i) = param;
    % Update to the lower bound of the perturbed parameter
    param_copies_l(i).(param_LSA_HPPC{i}) = param.(param_LSA_HPPC{i})*(1-pct);
    % Do not change code here (updating constrained parameters)
    %----------------------------------------------------------------------
    param_copies_l(i).ep_s=1-param_copies_l(i).ep;
    param_copies_l(i).en_s=1-param_copies_l(i).en;
    param_copies_l(i).ap=(3/param_copies_l(i).Rpp)*param_copies_l(i).ep_s;
    param_copies_l(i).an=(3/param_copies_l(i).Rpn)*param_copies_l(i).en_s;
    param_copies_l(i).curr_dens_vec = param_copies_l(i).curr_vec./param_copies_l(i).Acell;
    %----------------------------------------------------------------------

    %----------------------------------------------------------------------
    % Update param structures with the upper bound perturbed parameters
    %----------------------------------------------------------------------
    param_copies_u(i) = param;
    % Update to the upper bound of the perturbed parameter
    param_copies_u(i).(param_LSA_HPPC{i}) = param.(param_LSA_HPPC{i})*(1+pct);
    % Do not change code here (updating constrained parameters)
    %----------------------------------------------------------------------
    param_copies_u(i).ep_s=1-param_copies_u(i).ep;
    param_copies_u(i).en_s=1-param_copies_u(i).en;
    param_copies_u(i).ap=(3/param_copies_u(i).Rpp)*param_copies_u(i).ep_s;
    param_copies_u(i).an=(3/param_copies_u(i).Rpn)*param_copies_u(i).en_s;
    param_copies_u(i).curr_dens_vec = param_copies_u(i).curr_vec./param_copies_u(i).Acell;
    %----------------------------------------------------------------------
end

fprintf('Simulating with perturbed parameters...\n')
parfor (i = 1:length(param_LSA_HPPC), param.pso_workers)
% for i = 1:length(param_LSA_HPPC)
    %----------------------------------------------------------------------
    % Simulate using lower bound perturbed parameter
    %----------------------------------------------------------------------
    temp_params_l = param_copies_l(i);
    [t_l,V_l,SOCp_l,SOCn_l,~,~,~,~]=HPPC_sim(temp_params_l,temp_params_l.SOC_init,temp_params_l.t0,temp_params_l.tf_vec,temp_params_l.curr_dens_vec,temp_params_l.curr_segments);

    t_l_all{i} = t_l;
    V_l_all{i} = V_l;
    SOCp_l_all{i} = SOCp_l;
    SOCn_l_all{i} = SOCn_l;
    %----------------------------------------------------------------------
    % Simulate using upper bound perturbed parameter
    %----------------------------------------------------------------------
    temp_params_u = param_copies_u(i);
    [t_u,V_u,SOCp_u,SOCn_u,~,~,~,~]=HPPC_sim(temp_params_u,temp_params_u.SOC_init,temp_params_u.t0,temp_params_u.tf_vec,temp_params_u.curr_dens_vec,temp_params_u.curr_segments);

    t_u_all{i} = t_u;
    V_u_all{i} = V_u;
    SOCp_u_all{i} = SOCp_u;
    SOCn_u_all{i} = SOCn_u;
    %----------------------------------------------------------------------
    % Calculate voltage sensitivity vector
    %----------------------------------------------------------------------
    % Nominal parameter value
    theta = theta_nom(i); 
    % Choose minimum length to calculate sensitivity vectors
    len_V = min([length(V_u),length(V_nom),length(V_l)]);
    len_SOCp = min([length(SOCp_u),length(SOCp_nom),length(SOCp_l)]);
    len_SOCn = min([length(SOCn_u),length(SOCn_nom),length(SOCn_l)]);
    % Calculate normalized voltage sensitivity vector
    S_col_l_V = abs((V_l(1:len_V)-V_nom(1:len_V)))./abs(temp_params_l.(param_LSA_HPPC{i})-theta);
    S_col_u_V = abs((V_u(1:len_V)-V_nom(1:len_V)))./abs(temp_params_u.(param_LSA_HPPC{i})-theta);
    S_col_avg_V = (S_col_u_V+S_col_l_V)/2;
    S_col_norm_V = S_col_avg_V.*theta./V_nom(1:len_V);
    % Store calculated normalized voltage sensitivity vector in cell
    S_norm_V{i} = S_col_norm_V;
    % Calculate normalized SOCp sensitivity vector
    S_col_l_SOCp = abs((SOCp_l(1:len_SOCp)-SOCp_nom(1:len_SOCp)))./abs(temp_params_l.(param_LSA_HPPC{i})-theta);
    S_col_u_SOCp = abs((SOCp_u(1:len_SOCp)-SOCp_nom(1:len_SOCp)))./abs(temp_params_u.(param_LSA_HPPC{i})-theta);
    S_col_avg_SOCp = (S_col_l_SOCp+S_col_u_SOCp)/2;
    S_col_norm_SOCp = S_col_avg_SOCp.*theta./SOCp_nom(1:len_SOCp);
    % Store calculated normalized SOCp sensitivity vector in cell
    S_norm_SOCp{i} = S_col_norm_SOCp;
    % Calculate normalized SOCn sensitivity vector
    S_col_l_SOCn = abs((SOCn_l(1:len_SOCn)-SOCn_nom(1:len_SOCn)))./abs(temp_params_l.(param_LSA_HPPC{i})-theta);
    S_col_u_SOCn = abs((SOCn_u(1:len_SOCn)-SOCn_nom(1:len_SOCn)))./abs(temp_params_u.(param_LSA_HPPC{i})-theta);
    S_col_avg_SOCn = (S_col_l_SOCn+S_col_u_SOCn)/2;
    S_col_norm_SOCn = S_col_avg_SOCn.*theta./SOCn_nom(1:len_SOCn);
    % Store calculated normalized SOCn sensitivity vector in cell
    S_norm_SOCn{i} = S_col_norm_SOCn;
    % Store calculated normalized total sensitivity vector in cell
    S_norm_V_SOC{i} = [S_col_norm_V S_col_norm_SOCp S_col_norm_SOCn];
    fprintf(['Done perturbing ' param_LSA_HPPC{i} '...\n'])
end

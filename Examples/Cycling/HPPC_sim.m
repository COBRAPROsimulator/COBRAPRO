function [t_sim,V_sim,SOCp_sim,SOCn_sim,curr_dens_sim,...
    delta_phie_sim,delta_eta_sim,delta_OCV_sim,init_fail]=HPPC_sim(param,SOC_init,t0,tf_vec,currentDensity_vec,curr_segments)
%   HPPC_sim conducts HPPC simulation by conducting each CC segment in the 
%   HPPC profile and concatenating the results
%
%   Inputs:
%       param: parameter structure
%       SOC_init: initial SOC [-]
%
%   Output:
%       t_sim: simulated HPPC time vector [s]
%       V_sim: simulated voltage vector [V]
%       SOCp_sim: simulated positive electrode SOC vector [-]
%       SOCn_sim: simulated negative electrode SOC vector [-]
%       curr_dens_sim: simulated current density vector [A.m-2]
%       delta_phie_sim: simulated electrolyte potential difference between x=0 (CC|positive) and x=L (negative|CC) [V]
%       delta_eta_sim: simulated overpotential difference between x=0 (CC|PE) and x=L (NE|CC) [V]
%       delta_OCV_sim: simulated OCV difference between x=0 (CC|PE) and x=L (NE|CC) [V] [V]
%       init_fail: number indicating which CC segment failed to initialize (if init_fail = 0, all CC segments initialized successfully)

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

% Processing simulation settings
param.SOC_init = SOC_init;
param.t0 = t0;
param.tf_vec = tf_vec;
param.curr_dens_vec = currentDensity_vec;
param.curr_segments = curr_segments;

% Calculate total number of HPPC simulation data points taking into account user input param.deltaT_sim
HPPC_len = 0;
for i = 1:length(tf_vec)
    HPPC_len = length(0:param.deltaT_sim:tf_vec(i)) + HPPC_len;
end

% Preallocate variables for efficient memory
t_sim = zeros(1,HPPC_len);
V_sim = zeros(1,HPPC_len);
SOCp_sim = zeros(1,HPPC_len);
SOCn_sim = zeros(1,HPPC_len);
curr_dens_sim = zeros(1,HPPC_len);
delta_phie_sim = zeros(1,HPPC_len);
delta_eta_sim = zeros(1,HPPC_len);
delta_OCV_sim = zeros(1,HPPC_len);
phien_ind = param.Np+param.Ns+1:param.Np+param.Ns+param.Nn;

% Counter to keep track of which CC segment we are simulating in HPPC 
ctr = 1;
x_init = [];
%----------------------------------------------------------------------
% Start HPPC simulation
%----------------------------------------------------------------------
for j = 1:param.curr_segments
    % Extract current density and tf for this CC segment
    currentDensity = param.curr_dens_vec(j);
    tf = param.tf_vec(j);
    % Run CC segment simulation
    [output,param] = runModel(x_init,param,param.t0,tf,param.SOC_init,currentDensity);
    if param.sim_print == 1
        fprintf(['CC segment ',num2str(j),' done.\n'])
    end
    %----------------------------------------------------------------------
    % Successful initialization at the given CC segment
    %----------------------------------------------------------------------
    if ~isnan(output.t) 
        x_init = output.x_initials;
        V_sim(ctr:ctr+length(output.voltage)-1) = output.voltage;
        curr_dens_sim(ctr:ctr+length(output.curr_dens)-1) = output.curr_dens;
        SOCp_sim(ctr:ctr+length(output.soc_p_bulk)-1) = output.soc_p_bulk;
        SOCn_sim(ctr:ctr+length(output.soc_n_bulk)-1) = output.soc_n_bulk;
        phie_pcc = electrodeBC_hermite_interp(output.phie(:,1:param.Np),output.curr_dens,param,'p');
        phie_ncc = electrodeBC_hermite_interp(output.phie(:,phien_ind),output.curr_dens,param,'n');
        delta_phie_sim(ctr:ctr+length(output.voltage)-1) = phie_pcc - phie_ncc;
        eta_pcc = electrodeBC_hermite_interp(output.etap,output.curr_dens,param,'p');
        eta_ncc = electrodeBC_hermite_interp(output.etan,output.curr_dens,param,'n');
        delta_eta_sim(ctr:ctr+length(output.voltage)-1) = eta_pcc - eta_ncc;
        Upcc = electrodeBC_hermite_interp(output.Up,output.curr_dens,param,'p');
        Uncc = electrodeBC_hermite_interp(output.Un,output.curr_dens,param,'n');
        delta_OCV_sim(ctr:ctr+length(output.voltage)-1) = Upcc - Uncc;
        if j == 1 % For the first cycle
            t_sim(ctr:ctr+length(output.t)-1) = output.t;
            ctr = ctr + length(output.t);
        else % After the first cycle
            t_sim(ctr:ctr+length(output.t)-1) = max(t_sim) + 1 + output.t;
            ctr = ctr + length(output.t);
        end
        %----------------------------------------------------------------------
        % If any exit conditions are trigged, break out of HPPC simulation for loop
        %----------------------------------------------------------------------
        % Exit condition can be any of the following:
        % -> upper or lower voltage cut-offs
        % -> csp_max < csp_surf or csn_max < csn_surf (BV equation will be imaginary)
        % -> IDA solver failed during simulation for some other reason
        %----------------------------------------------------------------------
        if ~strcmp(output.exit_indicator,'NaN')
            % Since simulated ended preemptively, break out of HPPC simulation for loop
            break
        end
    %----------------------------------------------------------------------
    % Unsuccessful initialization at the given CC segment
    %----------------------------------------------------------------------
    else 
        % Since simulated ended preemptively, break out of HPPC simulation for loop
        break 
    end
end
% Indicate which CC segment failed if HPPC simulated failed due to initialization failure at a CC segment
init_fail = 0;
if isnan(output.t) 
    init_fail = j;
end

% If simulation ended preemptively, remove zeros and NaN 
nonzero_ind = t_sim~=0; 
nonzero_ind(1) = true;
t_sim = t_sim(nonzero_ind);
V_sim = V_sim(nonzero_ind);
SOCp_sim = SOCp_sim(nonzero_ind);
SOCn_sim = SOCn_sim(nonzero_ind);
delta_phie_sim = delta_phie_sim(nonzero_ind);
delta_eta_sim = delta_eta_sim(nonzero_ind);
delta_OCV_sim = delta_OCV_sim(nonzero_ind);
curr_dens_sim = curr_dens_sim(nonzero_ind);
V_sim=V_sim(~isnan(V_sim));
t_sim=t_sim(~isnan(V_sim));
SOCp_sim=SOCp_sim(~isnan(V_sim));
SOCn_sim=SOCn_sim(~isnan(V_sim));
curr_dens_sim=curr_dens_sim(~isnan(V_sim));
delta_phie_sim = delta_phie_sim(~isnan(delta_phie_sim));
delta_eta_sim = delta_eta_sim(~isnan(delta_eta_sim));
delta_OCV_sim = delta_OCV_sim(~isnan(delta_OCV_sim));

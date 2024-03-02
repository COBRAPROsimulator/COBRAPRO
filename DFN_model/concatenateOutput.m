function output_sim = concatenateOutput(t, x, param, output_sim)
%   concatenateOutput extracts the variables in x and stores in output_sim
%
%   Inputs:
%       t: scalar time [s]
%       x: differential and algebraic state vector
%       output_sim: output structure
%
%   Output:
%       output_sim: output structure

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

% Voltage calculation
phisp = x(param.phisp_ind);
phisn = x(param.phisn_ind);
switch param.VinterpMode
    case 1 % center of CV value 
        voltage = phisp(1) - phisn(end) - param.Rc*(-x(param.I_density_ind))*param.Acell;
    case 2 % linear interpolation
        phis_pcc = 1.5*phisp(1)-0.5*phisp(2);
        phis_ncc = 1.5*phisn(end)-0.5*phisn(end-1);
        voltage = phis_pcc - phis_ncc - param.Rc*(-x(param.I_density_ind))*param.Acell;
    case 3 % 3rd order hermite interpolation
        phis_pcc = electrodeBC_hermite_interp(phisp,x(param.I_density_ind),param,'p');
        phis_ncc = electrodeBC_hermite_interp(phisn,x(param.I_density_ind),param,'n');
        voltage = phis_pcc - phis_ncc - param.Rc*(-x(param.I_density_ind))*param.Acell;
end

% Extract variables from x and store in output_sim.
output_sim.t = [output_sim.t; t];
output_sim.x = [output_sim.x; x];
output_sim.voltage = [output_sim.voltage; voltage];
output_sim.ce = [output_sim.ce; x(param.ce_ind)*param.c0];
output_sim.csp = [output_sim.csp; x(param.csp_ind)*param.ctp];
output_sim.csn = [output_sim.csn; x(param.csn_ind)*param.ctn];
output_sim.phisp = [output_sim.phisp; phisp];
output_sim.phisn = [output_sim.phisn; phisn];
output_sim.phie = [output_sim.phie; x(param.phie_ind)];
output_sim.jp = [output_sim.jp; x(param.jp_ind)];
output_sim.jn = [output_sim.jn; x(param.jn_ind)];
output_sim.curr_dens = [output_sim.curr_dens; x(param.I_density_ind)];
% Extract surface, average, bulk particle concentration and electrode SOC 
if strcmp(param.cs_discret,'FVM')
    [csp_surf, csp_avg, csp_avg_bulk, soc_p_bulk,...
        csn_surf, csn_avg, csn_avg_bulk, soc_n_bulk] = ...
        calculateSolidConcentrationValues_FVM(x(param.csp_ind)*param.ctp,x(param.csn_ind)*param.ctn,x(param.jp_ind),x(param.jn_ind),param);
elseif strcmp(param.cs_discret,'FDM')
    [csp_surf, csp_avg, csp_avg_bulk, soc_p_bulk,...
        csn_surf, csn_avg, csn_avg_bulk, soc_n_bulk] = ...
        calculateSolidConcentrationValues_FDM(x(param.csp_ind)*param.ctp,x(param.csn_ind)*param.ctn,param);
end
output_sim.csp_surf = [output_sim.csp_surf; csp_surf];
output_sim.csp_avg = [output_sim.csp_avg; csp_avg];
output_sim.csp_avg_bulk = [output_sim.csp_avg_bulk; csp_avg_bulk];
output_sim.soc_p_bulk = [output_sim.soc_p_bulk; soc_p_bulk];
output_sim.csn_surf = [output_sim.csn_surf; csn_surf];
output_sim.csn_avg = [output_sim.csn_avg; csn_avg];
output_sim.csn_avg_bulk = [output_sim.csn_avg_bulk; csn_avg_bulk];
output_sim.soc_n_bulk = [output_sim.soc_n_bulk; soc_n_bulk];
output_sim.Up = [output_sim.Up; param.Up(csp_surf./param.ctp)];
output_sim.Un = [output_sim.Un; param.Un(csn_surf./param.ctn)];
phie = x(param.phie_ind);
output_sim.etap = [output_sim.etap; phisp - phie(1:param.Np) - param.Up(csp_surf./param.ctp)];
output_sim.etan = [output_sim.etan; phisn - phie(param.Np+param.Ns+1:param.Np+param.Ns+param.Nn) - param.Un(csn_surf./param.ctn)];
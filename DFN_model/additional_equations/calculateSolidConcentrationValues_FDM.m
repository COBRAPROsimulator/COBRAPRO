function [csp_surf_vec, csp_avg_vec, csp_avg_bulk, soc_p_bulk, csn_surf_vec, csn_avg_vec, csn_avg_bulk, soc_n_bulk] = calculateSolidConcentrationValues_FDM(csp, csn, param)
%   calculateSolidConcentrationValues calculates surface concentration,
%   average concentration, and bulk concentration, and SOC in the particles
%
%   Inputs:
%       csp: particle particle concentration in positive electrode [mol.m-3]
%       csn: particle particle concentration in negative electrode [mol.m-3]
%       param: Parameter structure
%
%   Output:
%       csp_surf_vec: Surface particle concentration in positive electrode [mol.m-3]
%       csp_avg_vec: Average particle concentration in positive electrode (along radial direction) [mol.m-3]
%       csp_avg_bulk: Bulk particle concentration in positive electrode (along x direction) [mol.m-3] 
%       soc_p_bulk: SOC in positive electrode [-]
%       csn_surf_vec: Surface particle concentration in negative electrode [mol.m-3]
%       csn_avg_vec: Average particle concentration in negative electrode (along radial direction) [mol.m-3]
%       csn_avg_bulk: Bulk particle concentration in negative electrode (along x direction) [mol.m-3]
%       soc_n_bulk: SOC in negative electrode [-]

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

% ------------------------------------------------------------------------
% Positive Electrode
% ------------------------------------------------------------------------

% Preallocate csp_surf for memory efficiency
csp_surf_vec = zeros(1,param.Np);
% Preallocate csp_avg for memory efficiency
csp_avg_vec = zeros(1,param.Np);
% positive average solid concentration along electrode direction
csp_avg_cv_integrated_total = zeros(1,param.Np);
% x-points for rectangle integration method for csp_avg_bulk
x_intval_p = linspace(0,param.lp,param.Np*2+1); 

for i=1:param.Np 
    r_first_index = (i-1)*(param.Nrp)+1; % radial position in x
    r_last_index = i*(param.Nrp); % radial position in x
    % cs,surf extract
    csp_surf_vec(i) = csp(r_last_index);
    % cs,avg calculation
    csp_center = 4/3*csp(r_first_index)-1/3*csp(r_first_index+1); % [Jangra 2019] 3-point forward-difference
    integrand = [csp_center,csp(r_first_index:r_last_index)].*(param.r_intval_p.^2);
    csp_avg_vec(i) = (3/param.Rpp^3)*trapz(param.r_intval_p,integrand,2); 
    % cs_avg_bulk
    csp_avg_integrand = ones(1,3).*csp_avg_vec(i); % populate cs value in 3 points for rectangle integration method
    csp_avg_cv_integrated = trapz(x_intval_p(i*2-1:i*2+1),csp_avg_integrand,2);
    csp_avg_cv_integrated_total(i) = csp_avg_cv_integrated;
end
% bulk SOC_p calculation
csp_avg_bulk_integrated = sum(csp_avg_cv_integrated_total,2); 
csp_avg_bulk = csp_avg_bulk_integrated/param.lp;
theta_p_bulk = csp_avg_bulk/param.ctp;
soc_p_bulk = (param.theta0_p - theta_p_bulk)/(param.theta0_p - param.theta100_p);

% make sure soc_bulk is not over 1 or below 0
if soc_p_bulk > 1
    soc_p_bulk = 1;
elseif soc_p_bulk < 0
    soc_p_bulk = 0;
end

% ------------------------------------------------------------------------
% Negative Electrode
% ------------------------------------------------------------------------

% Preallocate csn_surf for memory efficiency
csn_surf_vec = zeros(1,param.Nn);
% Preallocate csn_avg for memory efficiency
csn_avg_vec = zeros(1,param.Nn);
% Negative average solid concentration along electrode direction
csn_avg_cv_integrated_total = zeros(1,param.Nn);
% x-points for rectangle integration method for csp_avg_bulk
x_intval_n = linspace(0,param.ln,param.Nn*2+1); 

for i=1:param.Nn  
    r_first_index = (i-1)*(param.Nrn)+1; % radial position in x
    r_last_index = i*(param.Nrn); % radial position in x
    % cs,surf extract
    csn_surf_vec(i) = csn(r_last_index);
    % cs,avg calculation
    csn_center = 4/3*csn(r_first_index)-1/3*csn(r_first_index+1); % [Jangra 2019] 3-point forward-difference
    integrand = [csn_center,csn(r_first_index:r_last_index)].*(param.r_intval_n.^2);
    csn_avg_vec(i) = (3/param.Rpn^3)*trapz(param.r_intval_n,integrand,2); 
    % cs_avg_bulk
    csn_avg_integrand = ones(1,3).*csn_avg_vec(i); % populate cs value in 3 points for rectangle integration method
    csn_avg_cv_integrated = trapz(x_intval_n(i*2-1:i*2+1),csn_avg_integrand,2);
    csn_avg_cv_integrated_total(i) = csn_avg_cv_integrated;
end
% bulk SOC_n calculation
csn_avg_bulk_integrated = sum(csn_avg_cv_integrated_total,2); 
csn_avg_bulk = csn_avg_bulk_integrated/param.ln;
theta_n_bulk = csn_avg_bulk/param.ctn;
soc_n_bulk = (theta_n_bulk - param.theta0_n)/(param.theta100_n - param.theta0_n);

% make sure soc_bulk is not over 1 or below 0
if soc_n_bulk > 1
    soc_n_bulk = 1;
elseif soc_n_bulk < 0
    soc_n_bulk = 0;
end

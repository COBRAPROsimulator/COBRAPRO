function [csp_surf_vec, csp_avg_vec, csp_avg_bulk, soc_p_bulk, csn_surf_vec, csn_avg_vec, csn_avg_bulk, soc_n_bulk] = calculateSolidConcentrationValues_FVM(csp, csn, jp, jn, param)
%   calculateSolidConcentrationValues calculates surface concentration,
%   average concentration, and bulk concentration, and SOC in the particles
%
%   Inputs:
%       csp: particle particle concentration in positive electrode [mol.m-3]
%       csn: particle particle concentration in negative electrode [mol.m-3]
%       jp: positive pore-wall flux [mol.m-2.s-1]
%       jn: negative pore-wall flux [mol.m-2.s-1]
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
% Positive average solid concentration along electrode direction
csp_avg_cv_integrated_total = zeros(1,param.Np);
% x-points for rectangle integration method for csp_avg_bulk
x_intval_p = linspace(0,param.lp,param.Np*2+1); 
% r-points for rectangle integration method for csp_avg
r_intval_p = linspace(0,param.Rpp,param.Nrp*2+1); 

for i=1:param.Np 
    r_index_Nr = i*(param.Nrp); % last radial position in x
    % csp_surf calculation (3rd order Hermite interpolation)
    csp_outer_ghost = -jp(i)*param.rp/param.Dsp+csp(r_index_Nr); % csp value at r = Nrp+1
    csp_surf_vec(i) = cs_surf_hermite_interp(csp,r_index_Nr,csp_outer_ghost,param,'p');
    % Preallocate csp_avg_int for memory efficiency
    csp_avg_int_vec = zeros(1,param.Nrp);
    % csp_avg calculation
    for j = 1:param.Nrp
        r_index = (i-1)*param.Nrp+j;
        integrand = (ones(1,3).*csp(r_index)).*(r_intval_p(j*2-1:j*2+1).^2);
        % result from rectangular integration
        csp_avg_int_vec(j) = (3/param.Rpp^3)*trapz(r_intval_p(j*2-1:j*2+1),integrand,2);
    end
    csp_avg_vec(i) = sum(csp_avg_int_vec,2); 
    % cs_avg_bulk
    csp_avg_integrand = ones(1,3).*csp_avg_vec(i); % populate cs value in 3 points for rectangle integration method
    csp_avg_cv_integrated_total(i) = trapz(x_intval_p(i*2-1:i*2+1),csp_avg_integrand,2);
end
% bulk SOC_p calculation
csp_avg_bulk_integrated = sum(csp_avg_cv_integrated_total,2); 
csp_avg_bulk = csp_avg_bulk_integrated/param.lp;
theta_p_bulk = csp_avg_bulk/param.ctp;
soc_p_bulk = (param.theta0_p - theta_p_bulk)/(param.theta0_p - param.theta100_p);

% Make sure soc_bulk is not over 1 or below 0
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
x_intval_n = linspace(0,param.ln,param.Nn*2+1); % x-points for rectangle integration method for csp_avg_bulk
r_intval_n = linspace(0,param.Rpn,param.Nrn*2+1); % r-points for rectangle integration method for csn_avg

for i=1:param.Nn  
    r_index_Nr = i*(param.Nrn); % last radial position in x
    % csn_surf calculation (3rd order Hermite interpolation)
    csn_outer_ghost = -jn(i)*param.rn/param.Dsn+csn(r_index_Nr); % csn value at r = Nrn+1
    csn_surf_vec(i) = cs_surf_hermite_interp(csn,r_index_Nr,csn_outer_ghost,param,'n');
    % Preallocate csn_avg_int for memory efficiency
    csn_avg_int_vec = zeros(1,param.Nrn);
    % csn_avg calculation
    for j = 1:param.Nrn
        r_index = (i-1)*param.Nrn+j;
        integrand = (ones(1,3).*csn(r_index)).*(r_intval_n(j*2-1:j*2+1).^2);
        csn_avg_int_vec(j) = (3/param.Rpn^3)*trapz(r_intval_n(j*2-1:j*2+1),integrand,2);
    end
    csn_avg_vec(i) = sum(csn_avg_int_vec,2); 
    % cs_avg_bulk
    csn_avg_integrand = ones(1,3).*csn_avg_vec(i); % populate cs value in 3 points for rectangle integration method
    csn_avg_cv_integrated_total(i) = trapz(x_intval_n(i*2-1:i*2+1),csn_avg_integrand,2);
end
% bulk SOC_n calculation
csn_avg_bulk_integrated = sum(csn_avg_cv_integrated_total,2); 
csn_avg_bulk = csn_avg_bulk_integrated/param.ln;
theta_n_bulk = csn_avg_bulk/param.ctn;
soc_n_bulk = (theta_n_bulk - param.theta0_n)/(param.theta100_n - param.theta0_n);

% make sure SOC_n is not over 1 or below 0
if soc_n_bulk > 1
    soc_n_bulk = 1;
elseif soc_n_bulk < 0
    soc_n_bulk = 0;
end
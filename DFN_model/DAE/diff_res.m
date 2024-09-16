function [dx_tot, flag, new_data] = diff_res(t,x,xp,ida_user_data)
%   -> diff_res calculates the implicit ode system M(t,x,xp)=0 
%   if using the single-step method. M(t,x,xp)=0 consists of the perturbed
%   AEs (implicit odes) and the original ODEs with the switch function applied.
%
%   -> diff_res calculates the DAE in implicit form F(t,x,xp)=0
%   if using the SUNDIALS IDACalcIC method. F(t,x,xp)=0 comprises of the 
%   ODE and AE residual equations.
%
%   Inputs:
%       t: time
%       x: Vector consisting of algebraic and differential variables
%       xp: Vector consisting of derivative of algebraic and differential variables
%       ida_user_data: contains any functions or parameters defined by user
%
%   Output:
%       dx_tot: Vector where each entry corresponds to the DAE residual for each variable
%       flag: Required by IDA solver but not used in code
%       new_data: Required by IDA solver but not used in code

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

% This required to stop the simulation if the states are NaN. Without this, the simulation will get stuck.
if isnumeric(x)
    if isnan(x)
        return
    end
end

% Required by IDA solver but not used in this code
flag = 0;
new_data = [];

% Extract param structure from ida_user_data
param = ida_user_data.param;

import casadi.*
% Extract the node numbers for cleaner code
Np=param.Np;
Ns=param.Ns;
Nn=param.Nn;
Nrp=param.Nrp;
Nrn=param.Nrn;

%% Extract initialization method (Single-step or SUNDIALS IDACalcIC)
% -> If using single-step method, extract the perturbed algebraic equations
% (implicit ode equations) by calling the perturbed algebraic equation 
% function "res_perAlg" stored in ida_user_data.  
% -> If using IDACalcIC method, extract the algebraic residual equations 
% directly from function "alg_res.m".

% Extract initialization method from ida_user_data
init_Method = ida_user_data.init_Method;

% Algebraic variable indices within algebraic variable lineup
phisp_last_ind = param.phisp_length;
phisn_last_ind = phisp_last_ind + param.phisn_length;
phie_last_ind = phisn_last_ind + param.phie_length;
jp_last_ind = phie_last_ind + param.jp_length;
jn_last_ind = jp_last_ind + param.jn_length;

%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%% For single-step method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
if strcmp(init_Method,'SS')
    %----------------------------------------------------------------------
    % Calculate implicit odes in symbolic or numeric format (perturbed algebraic equations)
    %----------------------------------------------------------------------
    f_res_perAlg = ida_user_data.res_perAlg;
    res_perAlg = full(f_res_perAlg(x,xp));
    %----------------------------------------------------------------------
    % Extract implicit odes according to the variables.
    %----------------------------------------------------------------------
    % Extract positive solid phase potential.
    res_perAlg_phisp = res_perAlg(1:phisp_last_ind);
    % Extract negative solid phase potential.
    res_perAlg_phisn = res_perAlg(phisp_last_ind+1:phisn_last_ind);
    % Extract electrolyte potential.
    res_perAlg_phie = res_perAlg(phisn_last_ind+1:phie_last_ind);
    % Extract positive pore-wall flux.
    res_perAlg_jp = res_perAlg(phie_last_ind+1:jp_last_ind);
    % Extract positive pore-wall flux.
    res_perAlg_jn = res_perAlg(jp_last_ind+1:jn_last_ind);
    % Extract current density.
    if isnumeric(t) 
        param.I_density = param.currentDensity(t,param); 
    else
        param.I_density = 0; % dummy: since we are taking the Jacobian, doesn't matter what I_density is
    end
    % Perturb I_density to make into implicit ODE
    res_perAlg_curr_dens = param.mu*xp(param.I_density_ind)+x(param.I_density_ind)-param.I_density;   
    %----------------------------------------------------------------------
    % A single-step iteration-free initialization approach (Patented by Venkat Subramanian Group at UT Austin)- Allows for Only Acamedic purpose
    %----------------------------------------------------------------------
    ff=1/2*tanh(param.q*(t-param.initime))+1/2; 

%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%% For SUNDIALS IDACalcIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
elseif strcmp(init_Method,'IDACalcIC')
    %----------------------------------------------------------------------
    % Calculate algebraic residuals
    %----------------------------------------------------------------------
    x_tot = alg_res(x, param);
    %---------------------alg-------------------------------------------------
    % Extract AEs for each algebraic variable
    %----------------------------------------------------------------------
    % Extract positive solid phase potential
    res_perAlg_phisp = x_tot(1:phisp_last_ind);
    % Extract negative solid phase potential
    res_perAlg_phisn = x_tot(phisp_last_ind+1:phisn_last_ind);
    % Extract electrolyte potential
    res_perAlg_phie = x_tot(phisn_last_ind+1:phie_last_ind);
    % Extract positive pore-wall flux
    res_perAlg_jp = x_tot(phie_last_ind+1:jp_last_ind);
    % Extract positive pore-wall flux
    res_perAlg_jn = x_tot(jp_last_ind+1:jn_last_ind);
    % Extract current density
    if isnumeric(t) 
        param.I_density = param.currentDensity(t,param); 
    else
        param.I_density = 0; % dummy: since we are taking the Jacobian, doesn't matter what I_density is
    end
    res_perAlg_curr_dens = x(param.I_density_ind)-param.I_density;
    %----------------------------------------------------------------------
    % Switch function is set to 1 since not using single-step method
    %----------------------------------------------------------------------
    ff = 1;
end

%% Electrolyte concentration (-)

% Extract electrolyte concentration variables.
ce = x(param.ce_ind);
% Extract electrolyte concentration derivative variables.
dce = xp(param.ce_ind);
% Extract pore-wall flux variables (also used for csp and csn equations below).
jp = x(param.jp_ind);
jn = x(param.jn_ind);

% Calculate electrolyte diffusion coefficient at center of FVM, e.g., i = 1, 2, 3,..., N
De_eff_p = param.De_eff(ce(1:Np), param, 'p');
De_eff_s = param.De_eff(ce(Np+1:Np+Ns), param, 's');
De_eff_n = param.De_eff(ce(Np+Ns+1:end), param, 'n');

% Interpolate elecrolyte diffusion coefficient to find values at interface, e.g., i=1.5, 2.5, 3.5,..., N-0.5
[De_eff_p_inter, De_eff_s_inter, De_eff_n_inter] = interpolate_De_eff(De_eff_p,De_eff_s,De_eff_n,param);
De_eff_inter = [De_eff_p_inter; De_eff_s_inter; De_eff_n_inter];

% Calculate transference coefficient at center of FVM, e.g., i = 1, 2, 3,..., N
t1 = param.t1(ce, param);

% Preallocate the ce residual for memory efficiency.
if isnumeric(x)
    res_dce = zeros(length(ce),1);
else
    res_dce = SX.sym('res_ce',length(ce),1);
end

% ------------------------------------------------------------------------
% Positive Electrode
% ------------------------------------------------------------------------
hp_hs_median = (param.hp + param.hs)/2;

for i=1:Np            
    if i == 1 % Boundary equation (CC|positive)
        dce_i_flux_right = De_eff_inter(1)*(ce(2)-ce(1))/param.hp;
        dce_i_flux_left = 0; % d(ce(x=0))/dx = 0
    elseif i == Np % Interface equation (positive|separator)
        dce_i_flux_right = De_eff_inter(i)*(ce(i+1)-ce(i))/hp_hs_median;
        dce_i_flux_left = De_eff_inter(i-1)*(ce(i)-ce(i-1))/param.hp;
    else % Internal equation
        dce_i_flux_right = De_eff_inter(i)*(ce(i+1)-ce(i))/param.hp;
        dce_i_flux_left = De_eff_inter(i-1)*(ce(i)-ce(i-1))/param.hp;
    end
    dce_eq_rhs = (1/param.hp)*(dce_i_flux_right - dce_i_flux_left) + ...
                 param.ap*(1-t1(i))*jp(i)/param.c0;
    res_dce(i) = param.ep*dce(i) - ff*(dce_eq_rhs);
end  

% ------------------------------------------------------------------------
% Separator
% ------------------------------------------------------------------------
hs_hn_median = (param.hs + param.hn)/2;

for i=Np+1:Np+Ns 
    if i == Np+1 % Interface equation (positive|separator)
        dce_i_flux_right = De_eff_inter(i)*(ce(i+1)-ce(i))/param.hs;
        dce_i_flux_left = De_eff_inter(i-1)*(ce(i)-ce(i-1))/hp_hs_median;
    elseif i == Np+Ns % Interface equation (separator|negative)
        dce_i_flux_right = De_eff_inter(i)*(ce(i+1)-ce(i))/hs_hn_median;
        dce_i_flux_left = De_eff_inter(i-1)*(ce(i)-ce(i-1))/param.hs;
    else % Internal equation
        dce_i_flux_right = De_eff_inter(i)*(ce(i+1)-ce(i))/param.hs;
        dce_i_flux_left = De_eff_inter(i-1)*(ce(i)-ce(i-1))/param.hs;
    end
    dce_eq_rhs = (1/param.hs)*(dce_i_flux_right - dce_i_flux_left);
    res_dce(i) = param.es*dce(i) - ff*(dce_eq_rhs);
end

% ------------------------------------------------------------------------
% Negative Electrode
% ------------------------------------------------------------------------

for i=Np+Ns+1:Np+Ns+Nn
    if i == Np+Ns+1 % Interface equation (separator|negative)
        dce_i_flux_right = De_eff_inter(i)*(ce(i+1)-ce(i))/param.hn; % dce_sn_flux_right
        dce_i_flux_left = De_eff_inter(i-1)*(ce(i)-ce(i-1))/hs_hn_median; % dce_sn_flux_left
    elseif i == Np+Ns+Nn % Boundary equation (negative|CC)
        dce_i_flux_right = 0; 
        dce_i_flux_left = De_eff_inter(i-1)*(ce(i)-ce(i-1))/param.hn;
    else % Internal equation
        dce_i_flux_right = De_eff_inter(i)*(ce(i+1)-ce(i))/param.hn;
        dce_i_flux_left = De_eff_inter(i-1)*(ce(i)-ce(i-1))/param.hn;
    end
    flux_index = i - (Np+Ns); % Important!!
    dce_eq_rhs = (1/param.hn)*(dce_i_flux_right - dce_i_flux_left) + ...
                 param.an*(1-t1(i))*jn(flux_index)/param.c0;
    res_dce(i) = param.en*dce(i) - ff*(dce_eq_rhs);
end

%% Positive particle concentration (-)

% Extract positive solid concentration variables.
csp = x(param.csp_ind);
% Extract positive solid concentration derivative variables.
dcsp = xp(param.csp_ind);

% Preallocate the csp residual for memory efficiency.
if isnumeric(x)
    res_dcsp = zeros(length(csp),1);
else
    res_dcsp = SX.sym('res_dcsp',length(csp),1);
end

if strcmp(param.cs_discret,'FVM')
    % Governing equations
    for i=1:Np % index corresponding to x-position of particle  
        for j=1:Nrp % radial index 
            r_index = (i-1)*Nrp+j; % global r index at r = Nrp
            if j==Nrp % using BC at r=Rpp: dcs/dr(r=Rpp)=-jp/Dsp
                csp_outer_ghost = -jp(i)*param.rp/(param.Dsp*param.ctp)+csp(r_index); % csp value at r = Nrp+1
                res_dcsp(r_index) = dcsp(r_index) - ff*((3*param.Dsp/(param.rp*((j*param.rp)^3-((j-1)*param.rp)^3)))*((csp_outer_ghost-csp(r_index))*(j*param.rp)^2-(csp(r_index)-csp(r_index-1))*((j-1)*param.rp)^2));
            elseif j==1 % using BC at r=0: dcs/dr(r=0)=0
                res_dcsp(r_index) = dcsp(r_index) - ff*((3*param.Dsp/param.rp^2)*(csp(r_index+1)-csp(r_index))); 
            else % Internal equation FVM
                res_dcsp(r_index) = dcsp(r_index) - ff*((3*param.Dsp/(param.rp*((j*param.rp)^3-((j-1)*param.rp)^3)))*((csp(r_index+1)-csp(r_index))*(j*param.rp)^2-(csp(r_index)-csp(r_index-1))*((j-1)*param.rp)^2));
            end
        end
    end
elseif strcmp(param.cs_discret,'FDM')
    % Governing equations
    for i=1:Np % index corresponding to x-position of particle  
        for j=1:Nrp % radial index 
            r_index = (i-1)*Nrp+j; % index corresponding to r-position in csp
            if j==Nrp
                % [Allam] particle surface
                csp_outerGhost = csp(r_index-1)-2*param.rp*jp(i)/param.Dsp/param.ctp;
                res_dcsp(r_index) = dcsp(r_index) - ff*((param.Dsp/param.rp^2)*((j+1)/j*csp_outerGhost-2*csp(r_index)+(j-1)/j*csp(r_index-1)));
            elseif j==1
                res_dcsp(r_index) = dcsp(r_index) - ff*((param.Dsp/param.rp^2)*(2*csp(r_index+1)-2*csp(r_index)));
            else % Internal equation 
                res_dcsp(r_index) = dcsp(r_index) - ff*((param.Dsp/param.rp^2)*((j+1)/j*csp(r_index+1)-2*csp(r_index)+(j-1)/j*csp(r_index-1)));
            end
        end
    end
end

%% Negative particle concentration (-)

% Extract negative solid concentration variables.
csn = x(param.csn_ind);
% Extract negative solid concentration derivative variables.
dcsn = xp(param.csn_ind);

% Preallocate the csn residual for memory efficiency.
if isnumeric(x)
    res_dcsn = zeros(length(csn),1);
else
    res_dcsn = SX.sym('res_dcsn',length(csn),1);
end

if strcmp(param.cs_discret,'FVM')
    % Governing equations
    for i=1:Nn % index corresponding to x-position of particle  
        for j=1:Nrn % radial index 
            r_index = (i-1)*Nrn+j; % global r index at r = Nrn
            if j==Nrn % using BC at r=Rpn: dcs/dr(r=Rpn)=-jn/Dsn
                csn_outer_ghost = -jn(i)*param.rn/(param.Dsn*param.ctn)+csn(r_index); % csp value at r = Nrp+1
                res_dcsn(r_index) = dcsn(r_index) - ff*((3*param.Dsn/(param.rn*((j*param.rn)^3-((j-1)*param.rn)^3)))*((csn_outer_ghost-csn(r_index))*(j*param.rn)^2-(csn(r_index)-csn(r_index-1))*((j-1)*param.rn)^2));
            elseif j==1 % using BC at r=0: dcs/dr(r=0)=0
                res_dcsn(r_index) = dcsn(r_index) - ff*((3*param.Dsn/param.rn^2)*(csn(r_index+1)-csn(r_index))); 
            else % Internal equation FVM
                res_dcsn(r_index) = dcsn(r_index) - ff*((3*param.Dsn/(param.rn*((j*param.rn)^3-((j-1)*param.rn)^3)))*((csn(r_index+1)-csn(r_index))*(j*param.rn)^2-(csn(r_index)-csn(r_index-1))*((j-1)*param.rn)^2));
            end
        end
    end
elseif strcmp(param.cs_discret,'FDM')
    % Governing equations
    for i=1:Nn % index corresponding to x-position of particle       
        for j=1:Nrn % radial position 
            r_index = (i-1)*Nrn+j; % index corresponding to r-position in x
            if j==Nrn
                % [Allam] particle surface
                csn_outerGhost = csn(r_index-1)-2*param.rn*jn(i)/param.Dsn/param.ctn;
                res_dcsn(r_index) = dcsn(r_index) - ff*((param.Dsn/param.rn^2)*((j+1)/j*csn_outerGhost-2*csn(r_index)+(j-1)/j*csn(r_index-1))); 
            elseif j==1
                res_dcsn(r_index) = dcsn(r_index) - ff*((param.Dsn/param.rn^2)*(2*csn(r_index+1)-2*csn(r_index)));
            else % Internal equation 
                res_dcsn(r_index) = dcsn(r_index) - ff*((param.Dsn/param.rn^2)*((j+1)/j*csn(r_index+1)-2*csn(r_index)+(j-1)/j*csn(r_index-1)));
            end
        end  
    end
end

%% Build residual of differential equations.

% Assign residual of ODEs and AEs in correct sequence.
dx_tot = [res_dce;res_dcsp;res_dcsn;res_perAlg_phisp;res_perAlg_phisn;res_perAlg_phie;res_perAlg_jp;res_perAlg_jn;res_perAlg_curr_dens];

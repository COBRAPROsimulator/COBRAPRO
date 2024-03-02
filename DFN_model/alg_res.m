function x_tot = alg_res(x, param)
%   alg_res returns the algebraic residual equations, which consist
%   of the phis_p, phis_n, phie, jp, and jn governing equations.
%
%   The algebraic residual equations are given as:
%   g(t,x) = 0 where x includes both algebraic and differential variables.
%
%   Inputs:
%       x: Vector consisting of algebraic and differential variables
%       param: Parameter structure
%
%   Output:
%       x_tot: Vector where each entry corresponds to the algebraic equation residual for each algebraic variable
%           -> x_tot = g(t,x)

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

%% Extract variables. Use these variables to create the algebraic residual equations.
import casadi.*

% Extract electrolyte concentration variables.
ce = x(param.ce_ind);
% Extract positive solid concentration variables.
csp = x(param.csp_ind);
% Extract negative solid concentration variables.
csn = x(param.csn_ind);
% Extract solid phase potential variables.
phisp = x(param.phisp_ind);
phisn = x(param.phisn_ind);
% Extract liquid phase potential variables.
phie = x(param.phie_ind);
% Extract pore-wall flux variables.
jp = x(param.jp_ind);
jn = x(param.jn_ind);
% Extract current density variable.
I_density = x(param.I_density_ind);

% Extract the node numbers for cleaner code.
Np=param.Np;
Ns=param.Ns;
Nn=param.Nn;
Nrp=param.Nrp;
Nrn=param.Nrn;

%% Solid Potential (V)

% Preallocate the phisp residual for memory efficiency.
if isnumeric(x)
    res_phisp = zeros(length(phisp),1);
else
    res_phisp = SX.sym('res_phisp',length(phisp),1);
end

% Preallocate the phisn residual for memory efficiency.
if isnumeric(x)
    res_phisn = zeros(length(phisn),1);
else
    res_phisn = SX.sym('res_phisn',length(phisn),1);
end

% ------------------------------------------------------------------------
% Positive Electrode
% ------------------------------------------------------------------------
for i=1:Np                
    if i == 1 % Boundary equation (CC|positive)    
        phisp_flux_i_right = param.sigma_effp*(phisp(i+1)-phisp(i))/param.hp;
        phisp_flux_i_left = -I_density;
    elseif i == Np % Boundary equation (positive|separator)
        phisp_flux_i_right = 0; % phisp_flux_ps_right
        phisp_flux_i_left = param.sigma_effp*(phisp(i)-phisp(i-1))/param.hp;
    else % Internal equation (positive)
        phisp_flux_i_right = param.sigma_effp*(phisp(i+1)-phisp(i))/param.hp;
        phisp_flux_i_left = param.sigma_effp*(phisp(i)-phisp(i-1))/param.hp;
    end
    res_phisp(i) = phisp_flux_i_right - phisp_flux_i_left - param.hp*param.ap*param.F*jp(i);
end

% ------------------------------------------------------------------------
% Negative Electrode
% ------------------------------------------------------------------------
for i=1:Nn-1                
    if i == 1 % Boundary equation (separator|negative)
        phisn_flux_i_right = param.sigma_effn*(phisn(i+1)-phisn(i))/param.hn;
        phisn_flux_i_left = 0;
    else % Internal equation (negative)
        phisn_flux_i_right = param.sigma_effn*(phisn(i+1)-phisn(i))/param.hn;
        phisn_flux_i_left = param.sigma_effn*(phisn(i)-phisn(i-1))/param.hn;
    end
    res_phisn(i) = phisn_flux_i_right - phisn_flux_i_left - param.hn*param.an*param.F*jn(i);
end

% Boundary equation (negative|CC) 
i = Nn;
% Third order Hermite interpolation to write equation for phisn_cc
% x_Nn = param.hn*(Nn-0.5); % x distance from i = 0.5 to i = Nn
% x_Nn_1 = param.hn*(Nn-1.5); % x distance from i = 0.5 to i = Nn-1
% x1 = (1+2*(param.ln-x_Nn_1)/(x_Nn-x_Nn_1))*phisn(i-1)+(param.ln-x_Nn_1)*((phisn(i)-phisn(i-2))/(2*param.hn));
% x2 = ((param.ln-x_Nn)/(x_Nn_1-x_Nn))^2;
% x3 = (1+2*(x_Nn-param.ln)/(x_Nn-x_Nn_1))*phisn(i)+(param.ln-x_Nn)*(-I_density/param.sigma_effn);
% x4 = ((param.ln-x_Nn_1)/(x_Nn-x_Nn_1))^2;
% phisn_end = x1*x2 + x3*x4;
phisn_end = electrodeBC_hermite_interp(phisn',I_density,param,'n');
res_phisn(i) = phisn_end;

%% Liquid phase potential (V)

% Preallocate the phie residual for memory efficiency.
if isnumeric(x)
    res_phie = zeros(length(phie),1);
else
    res_phie = SX.sym('res_phie',length(phie),1);
end

% ------------------------------------------------------------------------
% Calculate Kappa_eff, t1, lambda at center and interface of CVs (ce dependent coefficients)
% ------------------------------------------------------------------------
% Calculate electrolyte conductivity coefficient at center of FVM, e.g., i = 1, 2, 3,..., N
Kappa_eff_p = param.Kappa_eff(ce(1:Np), param, 'p');
Kappa_eff_s = param.Kappa_eff(ce(Np+1:Np+Ns), param, 's');
Kappa_eff_n = param.Kappa_eff(ce(Np+Ns+1:end), param, 'n');

% Interpolate elecrolyte diffusion coefficient to find values at interface, e.g., i=1.5, 2.5, 3.5,..., N-0.5
[Kappa_eff_p_inter, Kappa_eff_s_inter, Kappa_eff_n_inter] = interpolate_Kappa_eff(Kappa_eff_p,Kappa_eff_s,Kappa_eff_n,param);
Kappa_eff_inter = [Kappa_eff_p_inter; Kappa_eff_s_inter; Kappa_eff_n_inter];

% Calculate transference coefficient at center of FVM, e.g., i = 1, 2, 3,..., N
t1 = param.t1(ce, param);

% Interpolate transference coefficient to find values at interface, e.g., i=1.5, 2.5, 3.5,..., N-0.5
% i.e., i=1.5, 2.5, 3.5,...
[t1_p_inter, t1_s_inter, t1_n_inter] = interpolate_t1(t1,param);
t1_inter = [t1_p_inter; t1_s_inter; t1_n_inter];

% Calculate activity term at center of FVM, e.g., i = 1, 2, 3,..., N
lambda = param.lambda(ce, param);

% Interpolate activity term to find values at interface, e.g., i=1.5, 2.5, 3.5,..., N-0.5
[lambda_p_inter, lambda_s_inter, lambda_n_inter] = interpolate_lambda(lambda,param);
lambda_inter = [lambda_p_inter; lambda_s_inter; lambda_n_inter];

% Interpolate electrolyte concentration to find values at interface, e.g., i=1.5, 2.5, 3.5,..., N-0.5
[ce_p_inter, ce_s_inter, ce_n_inter] = interpolate_ce(ce,param);
ce_inter = [ce_p_inter; ce_s_inter; ce_n_inter];

% Define gamma parameter (constant value since we have isothermal conditions)
gamma = 2*param.T*param.R/param.F;

% ------------------------------------------------------------------------
% Positive Electrode
% ------------------------------------------------------------------------
hp_hs_median = (param.hp + param.hs)/2;

for i=1:Np         
    if i == 1 % Boundary equation (CC|positive)
        phie_flux_i_right = Kappa_eff_inter(i)*(phie(i+1)-phie(i))/param.hp;
        phie_flux_i_left = 0; % d(phie(x=0)/dx = 0
        log_ce_flux_i_right = Kappa_eff_inter(i)*gamma*(1-t1_inter(i))*lambda_inter(i)*(ce(i+1)-ce(i))/ce_inter(i)/param.hp;
        log_ce_flux_i_left = 0; % d(ce(x=0)/dx = 0
    elseif i == Np % Interface equation (positive|separator)
        phie_flux_i_right = Kappa_eff_inter(i)*(phie(i+1)-phie(i))/hp_hs_median; % phie_flux_ps_right
        phie_flux_i_left = Kappa_eff_inter(i-1)*(phie(i)-phie(i-1))/param.hp;
        log_ce_flux_i_right = Kappa_eff_inter(i)*gamma*(1-t1_inter(i))*lambda_inter(i)*(ce(i+1)-ce(i))/ce_inter(i)/hp_hs_median; % log_ce_flux_ps_right
        log_ce_flux_i_left = Kappa_eff_inter(i-1)*gamma*(1-t1_inter(i-1))*lambda_inter(i-1)*(ce(i)-ce(i-1))/ce_inter(i-1)/param.hp;
    else
        % Internal equation (positive)
        phie_flux_i_right = Kappa_eff_inter(i)*(phie(i+1)-phie(i))/param.hp;
        phie_flux_i_left = Kappa_eff_inter(i-1)*(phie(i)-phie(i-1))/param.hp;
        log_ce_flux_i_right = Kappa_eff_inter(i)*gamma*(1-t1_inter(i))*lambda_inter(i)*(ce(i+1)-ce(i))/ce_inter(i)/param.hp;
        log_ce_flux_i_left = Kappa_eff_inter(i-1)*gamma*(1-t1_inter(i-1))*lambda_inter(i-1)*(ce(i)-ce(i-1))/ce_inter(i-1)/param.hp;
    end
    res_phie(i) = -(phie_flux_i_right - phie_flux_i_left) + ...
               (log_ce_flux_i_right - log_ce_flux_i_left) - param.hp*param.ap*param.F*jp(i);
end

% ------------------------------------------------------------------------
% Separator
% ------------------------------------------------------------------------
hs_hn_median = (param.hs + param.hn)/2;

for i=Np+1:Np+Ns  
    if i == Np+1 % Interface equation (positive|separator)
        phie_flux_i_right = Kappa_eff_inter(i)*(phie(i+1)-phie(i))/param.hs;
        phie_flux_i_left = Kappa_eff_inter(i-1)*(phie(i)-phie(i-1))/hp_hs_median; % phie_flux_ps_left
        log_ce_flux_i_right = Kappa_eff_inter(i)*gamma*(1-t1_inter(i))*lambda_inter(i)*(ce(i+1)-ce(i))/ce_inter(i)/param.hs;
        log_ce_flux_i_left = Kappa_eff_inter(i-1)*gamma*(1-t1_inter(i-1))*lambda_inter(i-1)*(ce(i)-ce(i-1))/ce_inter(i-1)/hp_hs_median; % log_ce_flux_ps_left
    elseif i == Np+Ns % Interface equation (separator|negative)
        phie_flux_i_right = Kappa_eff_inter(i)*(phie(i+1)-phie(i))/hs_hn_median; % phie_flux_sn_right
        phie_flux_i_left = Kappa_eff_inter(i-1)*(phie(i)-phie(i-1))/param.hs;
        log_ce_flux_i_right = Kappa_eff_inter(i)*gamma*(1-t1_inter(i))*lambda_inter(i)*(ce(i+1)-ce(i))/ce_inter(i)/hs_hn_median; % log_ce_flux_sn_right
        log_ce_flux_i_left = Kappa_eff_inter(i-1)*gamma*(1-t1_inter(i-1))*lambda_inter(i-1)*(ce(i)-ce(i-1))/ce_inter(i-1)/param.hs;
    else
        % Internal equation (separator)
        phie_flux_i_right = Kappa_eff_inter(i)*(phie(i+1)-phie(i))/param.hs;
        phie_flux_i_left = Kappa_eff_inter(i-1)*(phie(i)-phie(i-1))/param.hs;
        log_ce_flux_i_right = Kappa_eff_inter(i)*gamma*(1-t1_inter(i))*lambda_inter(i)*(ce(i+1)-ce(i))/ce_inter(i)/param.hs;
        log_ce_flux_i_left = Kappa_eff_inter(i-1)*gamma*(1-t1_inter(i-1))*lambda_inter(i-1)*(ce(i)-ce(i-1))/ce_inter(i-1)/param.hs;
    end
    res_phie(i) = -(phie_flux_i_right - phie_flux_i_left) + ...
               (log_ce_flux_i_right - log_ce_flux_i_left);
end

% ------------------------------------------------------------------------
% Negative Electrode
% ------------------------------------------------------------------------

for i=Np+Ns+1:Np+Ns+Nn 
    if i == Np+Ns+1 % Interface equation (separator|negative)
        phie_flux_i_right = Kappa_eff_inter(i)*(phie(i+1)-phie(i))/param.hn;
        phie_flux_i_left = Kappa_eff_inter(i-1)*(phie(i)-phie(i-1))/hs_hn_median; % phie_flux_sn_left
        log_ce_flux_i_right = Kappa_eff_inter(i)*gamma*(1-t1_inter(i))*lambda_inter(i)*(ce(i+1)-ce(i))/ce_inter(i)/param.hn;
        log_ce_flux_i_left = Kappa_eff_inter(i-1)*gamma*(1-t1_inter(i-1))*lambda_inter(i-1)*(ce(i)-ce(i-1))/ce_inter(i-1)/hs_hn_median; % log_ce_flux_sn_left
    elseif i == Np+Ns+Nn % Boundary equation (negative|CC)
        phie_flux_i_right = 0; % d(phie(x=0)/dx = 0
        phie_flux_i_left = Kappa_eff_inter(i-1)*(phie(i)-phie(i-1))/param.hn;
        log_ce_flux_i_right = 0; % d(ce(x=0)/dx = 0
        log_ce_flux_i_left = Kappa_eff_inter(i-1)*gamma*(1-t1_inter(i-1))*lambda_inter(i-1)*(ce(i)-ce(i-1))/ce_inter(i-1)/param.hn;
    else
        % Internal equation
        phie_flux_i_right = Kappa_eff_inter(i)*(phie(i+1)-phie(i))/param.hn;
        phie_flux_i_left = Kappa_eff_inter(i-1)*(phie(i)-phie(i-1))/param.hn;
        log_ce_flux_i_right = Kappa_eff_inter(i)*gamma*(1-t1_inter(i))*lambda_inter(i)*(ce(i+1)-ce(i))/ce_inter(i)/param.hn;
        log_ce_flux_i_left = Kappa_eff_inter(i-1)*gamma*(1-t1_inter(i-1))*lambda_inter(i-1)*(ce(i)-ce(i-1))/ce_inter(i-1)/param.hn;
    end
    flux_index = i - (Np+Ns); % Important!! 
    res_phie(i) = -(phie_flux_i_right - phie_flux_i_left) + ...
               (log_ce_flux_i_right - log_ce_flux_i_left) - param.hn*param.an*param.F*jn(flux_index);
end

%% Pore-wall flux (mol.m-2.s-1)

% Preallocate the jp residual for memory efficiency.
if isnumeric(x)
    res_jp = zeros(length(jp),1);
else
    res_jp = SX.sym('res_jp',length(jp),1);
end

% Preallocate the jn residual for memory efficiency.
if isnumeric(x)
    res_jn = zeros(length(jn),1);
else
    res_jn = SX.sym('res_jn',length(jn),1);
end

% ------------------------------------------------------------------------
% Positive Electrode
% ------------------------------------------------------------------------
for i=1:Np % x-direction index referring to particle locations
    if strcmp(param.cs_discret,'FVM')
        r_index_Nrp = i*Nrp; % global r index at r = Nrp
        cs_outer_ghost = -jp(i)*param.rp/(param.Dsp*param.ctp)+csp(r_index_Nrp); % csp value at r = Nrp+1
        csp_surf = cs_surf_hermite_interp(csp,r_index_Nrp,cs_outer_ghost,param,'p');
        theta=csp_surf; % csp is already normalized by csp_max
    elseif strcmp(param.cs_discret,'FDM')
        surf_index = i*Nrp; % particle surface index
        theta=csp(surf_index); % csp is already normalized by csp_max
    end
    Up=param.Up(theta);
    jp0 = 2*param.kp*((ce(i)*param.c0).^0.5)*((param.ctp-theta*param.ctp).^0.5)*((theta*param.ctp).^0.5); % exchange current density
    res_jp(i) = jp(i) - jp0*sinh(0.5*param.F/param.R/param.T*(phisp(i)-phie(i)-Up));
end

% Exit function since surface particle concentration is greater than maximum concentration
if isnumeric(x)
    if param.ctp-theta*param.ctp < 0
        return
    end
end

% ------------------------------------------------------------------------
% Negative Electrode
% ------------------------------------------------------------------------
 
for i=1:Nn % local index referring to particle locations
    if strcmp(param.cs_discret,'FVM')
        r_index_Nrn = i*Nrn; % global r index at r = Nrn
        cs_outer_ghost = -jn(i)*param.rn/(param.Dsn*param.ctn)+csn(r_index_Nrn); % csp value at r = Nrn+1
        csn_surf = cs_surf_hermite_interp(csn,r_index_Nrn,cs_outer_ghost,param,'n');
        theta=csn_surf; % csn is already normalized by csn_max
    elseif strcmp(param.cs_discret,'FDM')
        surf_index = i*Nrn; % particle surface index
        theta=csn(surf_index); % csn is already normalized by csp_max
    end
    Un=param.Un(theta);
    electrode_index = Np+Ns+i;
    jn0 = 2*param.kn*((ce(electrode_index)*param.c0).^0.5)*((param.ctn-theta*param.ctn).^0.5)*((theta*param.ctn).^0.5); % exchange current density
    res_jn(i) = jn(i) - jn0*sinh(0.5*param.F/param.R/param.T*(phisn(i)-phie(electrode_index)-Un));
end

% Exit function since surface particle concentration is greater than maximum concentration
if isnumeric(x)
    if param.ctn-theta*param.ctn < 0
        return
    end
end

%% Build algebraic residual

% Concatenate rsidual of algebraic equations in the correct sequence
x_tot = [res_phisp; res_phisn; res_phie; res_jp; res_jn];
function [ce_p_inter, ce_s_inter, ce_n_inter] = interpolate_ce(ce,param)
%   interpolate_ce interpolates the value of electrolyte concentration 
%   at the CV interface (i.e., edges) using harmonic mean
%
%   Inputs:
%       ce: normalized electrolyte concentration [-]
%       param: Parameter structure
%
%   Output:
%       ce_p_inter: normalized electrolyte concentration at the positive electrode CV edges, s.t., i=1.5, 2.5, 3.5,..., Np+0.5 [-]
%       ce_s_inter: normalizedelectrolyte concentration at the separator CV edges, s.t., i=Np+1.5, 2.5, 3.5,..., Np+Ns+0.5 [-]
%       ce_n_inter: normalized electrolyte concentration at the negative electrode CV edges, s.t., i= Np+Ns+1.5, 2.5, 3.5,..., Np+Ns+Nn-0.5 [-]

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

%% Interpolation within the positive electrode
beta_ce_p = 0.5;
ce_mean_p = (ce(1:param.Np-1).*ce(2:param.Np)) ./ (beta_ce_p*ce(2:param.Np) + (1-beta_ce_p)*ce(1:param.Np-1));

% Interpolation on the interface between separator and positive electrode
beta_ce_ps = param.hp / (param.hp + param.hs);
ce_mean_ps = (ce(param.Np)*ce(param.Np+1)) / (beta_ce_ps*ce(param.Np+1) + (1-beta_ce_ps)*ce(param.Np));

% ce values at CV edges in positive electrode
ce_p_inter = [ce_mean_p; ce_mean_ps];

%% Interpolation within the separator
beta_ce_s = 0.5;
ce_mean_s = ce(param.Np+1:param.Np+param.Ns-1).*ce(param.Np+2:param.Np+param.Ns)./ (beta_ce_s*ce(param.Np+2:param.Np+param.Ns) + (1-beta_ce_s)*ce(param.Np+1:param.Np+param.Ns-1));

% Interpolation on the interface between separator and negative electrode
beta_ce_sn = param.hs / (param.hn + param.hs);
ce_mean_sn = ce(param.Np+param.Ns)*ce(param.Np+param.Ns+1) / (beta_ce_sn*ce(param.Np+param.Ns+1) + (1-beta_ce_sn)*ce(param.Np+param.Ns));

% ce values at CV edges in separator 
ce_s_inter = [ce_mean_s; ce_mean_sn];

%% Interpolation within the negative electrode
beta_ce_n = 0.5;
% ce values at CV edges in negative electrode
ce_n_inter = ce(param.Np+param.Ns+1:end-1).*ce(param.Np+param.Ns+2:end)./ (beta_ce_n*ce(param.Np+param.Ns+2:end) + (1-beta_ce_n)*ce(param.Np+param.Ns+1:end-1));

end
function [t1_p_inter, t1_s_inter, t1_n_inter] = interpolate_t1(t1,param)
%   interpolate_t1 interpolates the transference number at the CV interface 
%   (i.e., edges) using harmonic mean
%
%   Inputs:
%       t1: transference number [-]
%       param: Parameter structure
%
%   Output:
%       t1_p_inter: transference number at the positive electrode CV edges, s.t., i=1.5, 2.5, 3.5,..., Np+0.5 [-]
%       t1_s_inter: transference number at the separator CV edges, s.t., i=Np+1.5, 2.5, 3.5,..., Np+Ns+0.5 [-]
%       t1_n_inter: transference number at the negative electrode CV edges, s.t., i= Np+Ns+1.5, 2.5, 3.5,..., Np+Ns+Nn-0.5 [-]

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

%% Interpolation in positive electrode

% Interpolation within the positive electrode
beta_t1_p = 0.5;
t1_mean_p = (t1(1:param.Np-1).*t1(2:param.Np)) ./ (beta_t1_p*t1(2:param.Np) + (1-beta_t1_p)*t1(1:param.Np-1));

% Interpolation on the interface between separator and positive electrode
beta_t1_ps = param.hp / (param.hp + param.hs);
t1_mean_ps = (t1(param.Np)*t1(param.Np+1)) / (beta_t1_ps*t1(param.Np+1) + (1-beta_t1_ps)*t1(param.Np));

% t1_p at CV edges in positive electrode
t1_p_inter = [t1_mean_p; t1_mean_ps];

%% Interpolation in separator

% Interpolation within the separator
beta_t1_s = 0.5;
t1_mean_s = t1(param.Np+1:param.Np+param.Ns-1).*t1(param.Np+2:param.Np+param.Ns)./ (beta_t1_s*t1(param.Np+2:param.Np+param.Ns) + (1-beta_t1_s)*t1(param.Np+1:param.Np+param.Ns-1));

% Interpolation on the interface between negative electrode and separator
beta_t1_sn = param.hs / (param.hs + param.hn);
t1_mean_sn = (t1(param.Ns)*t1(param.Ns+1)) / (beta_t1_sn*t1(param.Ns+1) + (1-beta_t1_sn)*t1(param.Ns));

% t1_s at CV edges in separator
t1_s_inter = [t1_mean_s; t1_mean_sn];

%% Interpolation in negative electrode

% Interpolation within the negative electrode
% t1_n at CV edges in negative electrode
beta_t1_n = 0.5;
t1_n_inter = t1(param.Ns+1:param.Ns+param.Nn-1).*t1(param.Ns+2:param.Ns+param.Nn)./ (beta_t1_n*t1(param.Ns+2:param.Ns+param.Nn) + (1-beta_t1_n)*t1(param.Ns+1:param.Ns+param.Nn-1));

end

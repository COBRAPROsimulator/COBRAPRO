function [lambda_p_inter, lambda_s_inter, lambda_n_inter] = interpolate_lambda(lambda,param)
%   interpolate_lambda interpolates the value of the activity coefficient 
%   at the CV interface (i.e., edges) using harmonic mean
%
%   Inputs:
%       lambda: activity coefficient [-]
%       param: Parameter structure
%
%   Output:
%       lambda_p_inter: activity coefficient at the positive electrode CV edges, s.t., i=1.5, 2.5, 3.5,..., Np+0.5 [-]
%       lambda_s_inter: activity coefficient at the separator CV edges, s.t., i=Np+1.5, 2.5, 3.5,..., Np+Ns+0.5 [-]
%       lambda_n_inter: activity coefficient at the negative electrode CV edges, s.t., i= Np+Ns+1.5, 2.5, 3.5,..., Np+Ns+Nn-0.5 [-]

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
beta_lambda_p = 0.5;
lambda_mean_p = (lambda(1:param.Np-1).*lambda(2:param.Np)) ./ (beta_lambda_p*lambda(2:param.Np) + (1-beta_lambda_p)*lambda(1:param.Np-1));

% Interpolation on the interface between separator and positive electrode
beta_lambda_ps = param.hp / (param.hp + param.hs);
lambda_mean_ps = (lambda(param.Np)*lambda(param.Np+1)) / (beta_lambda_ps*lambda(param.Np+1) + (1-beta_lambda_ps)*lambda(param.Np));

% lambda_p at CV edges in positive electrode
lambda_p_inter = [lambda_mean_p; lambda_mean_ps];

%% Interpolation in separator

% Interpolation within the separator
beta_lambda_s = 0.5;
lambda_mean_s = lambda(param.Np+1:param.Np+param.Ns-1).*lambda(param.Np+2:param.Np+param.Ns)./ (beta_lambda_s*lambda(param.Np+2:param.Np+param.Ns) + (1-beta_lambda_s)*lambda(param.Np+1:param.Np+param.Ns-1));

% Interpolation on the interface between negative electrode and separator
beta_lambda_sn = param.hs / (param.hs + param.hn);
lambda_mean_sn = (lambda(param.Ns)*lambda(param.Ns+1)) / (beta_lambda_sn*lambda(param.Ns+1) + (1-beta_lambda_sn)*lambda(param.Ns));

% lambda_s at CV edges in separator
lambda_s_inter = [lambda_mean_s; lambda_mean_sn];

%% Interpolation in negative electrode

% Interpolation within the negative electrode
beta_lambda_n = 0.5;
lambda_mean_n = lambda(param.Ns+1:param.Ns+param.Nn-1).*lambda(param.Ns+2:param.Ns+param.Nn)./ (beta_lambda_n*lambda(param.Ns+2:param.Ns+param.Nn) + (1-beta_lambda_n)*lambda(param.Ns+1:param.Ns+param.Nn-1));

% lambda_n at CV edges in negative electrode
lambda_n_inter = lambda_mean_n;

end
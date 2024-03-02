function [Kappa_eff_p_inter, Kappa_eff_s_inter, Kappa_eff_n_inter] = interpolate_Kappa_eff(Kappa_eff_p,Kappa_eff_s,Kappa_eff_n,param)
%   interpolate_Kappa_eff interpolates the value of the effective electrolyte
%   conductivity coefficient at the CV interface (i.e., edges) using harmonic mean
%
%   Inputs:
%       Kappa_eff_p: effective electrolyte diffusion in positive electrode [S.m-1]
%       Kappa_eff_s: effective electrolyte diffusion in separator [S.m-1]
%       Kappa_eff_n: effective electrolyte diffusion in negative electrode [S.m-1]
%       param: Parameter structure
%
%   Output:
%       Kappa_eff_p_inter: effective electrolyte diffusion at the positive electrode CV edges, s.t., i=1.5, 2.5, 3.5,..., Np+0.5 [S.m-1]
%       Kappa_eff_s_inter: effective electrolyte diffusion at the separator CV edges, s.t., i=Np+1.5, 2.5, 3.5,..., Np+Ns+0.5 [S.m-1]
%       Kappa_eff_n_inter: effective electrolyte diffusion at the negative electrode CV edges, s.t., i= Np+Ns+1.5, 2.5, 3.5,..., Np+Ns+Nn-0.5 [S.m-1]

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

%% Positive electrode mean conductivity

% Beta coefficient within the positive electrode
beta_p = 0.5;
% Harmonic mean for the coefficients within the positive electrode
Kappa_eff_p_inter = (Kappa_eff_p(1:end-1).*Kappa_eff_p(2:end)) ./ (beta_p*Kappa_eff_p(2:end)+(1-beta_p)*Kappa_eff_p(1:end-1));
% Beta coefficient for the interface with the separator
beta_ps = param.hp / (param.hp + param.hs);
% Harmonic mean for the coefficients at the interface with the separator
Kappa_eff_inter_ps = (Kappa_eff_p(end)*Kappa_eff_s(1)) / (beta_ps*Kappa_eff_s(1) + (1-beta_ps)*Kappa_eff_p(end));
% Kappa_eff_p at CV edges in positive electrode
Kappa_eff_p_inter = [Kappa_eff_p_inter;Kappa_eff_inter_ps];

%% Separator mean conductivity

% Beta coefficient within the separator
beta_s = 0.5;
% Harmonic mean for the coefficients within the separator
Kappa_eff_s_inter = (Kappa_eff_s(1:end-1).*Kappa_eff_s(2:end)) ./ (beta_s*Kappa_eff_s(2:end)+(1-beta_s)*Kappa_eff_s(1:end-1));
% Beta coefficient for the interface with the negative electrode
betaD_sn = param.hs / (param.hs + param.hn);
% Harmonic mean for the coefficients at the interface with the negative electrode
Kappa_eff_inter_sn = (Kappa_eff_s(end)*Kappa_eff_n(1)) / (betaD_sn*Kappa_eff_n(1) + (1-betaD_sn)*Kappa_eff_s(end));
% Kappa_eff_s at CV edges in separeator
Kappa_eff_s_inter = [Kappa_eff_s_inter;Kappa_eff_inter_sn];

%% Negative electrode mean conductivity

% Beta coefficient within the negative electrode
betaD_n = 0.5;
% Kappa_eff_n at CV edges in negative electrode
Kappa_eff_n_inter = Kappa_eff_n(1:end-1).*Kappa_eff_n(2:end)./(betaD_n*Kappa_eff_n(2:end)+(1-betaD_n)*Kappa_eff_n(1:end-1));

end
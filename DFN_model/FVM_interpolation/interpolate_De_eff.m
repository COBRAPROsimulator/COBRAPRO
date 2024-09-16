function [De_eff_p_inter, De_eff_s_inter, De_eff_n_inter] = interpolate_De_eff(De_eff_p,De_eff_s,De_eff_n,param)
%   interpolate_De_eff interpolates the value of the effective electrolyte
%   diffusion coefficient at the CV interface (i.e., edges) using harmonic mean
%
%   Inputs:
%       De_eff_p: effective electrolyte diffusion in positive electrode [m2.s-1]
%       De_eff_s: effective electrolyte diffusion in separator [m2.s-1]
%       De_eff_n: effective electrolyte diffusion in negative electrode [m2.s-1]
%       param: Parameter structure
%
%   Output:
%       De_eff_p_inter: effective electrolyte diffusion at the positive electrode CV edges, s.t., i=1.5, 2.5, 3.5,..., Np+0.5 [m2.s-1]
%       De_eff_s_inter: effective electrolyte diffusion at the separator CV edges, s.t., i=Np+1.5, 2.5, 3.5,..., Np+Ns+0.5 [m2.s-1]
%       De_eff_n_inter: effective electrolyte diffusion at the negative electrode CV edges, s.t., i= Np+Ns+1.5, 2.5, 3.5,..., Np+Ns+Nn-0.5 [m2.s-1]

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

%% Diffusion coefficients for the positive electrode

% Beta coefficient within the positive electrode
betaD_p = 0.5;
% Harmonic mean for the coefficients within the positive electrode
De_eff_p_inter = (De_eff_p(1:end-1).*De_eff_p(2:end)) ./ (betaD_p*De_eff_p(2:end)+(1-betaD_p)*De_eff_p(1:end-1));
% Beta coefficient for the interface with the separator
betaD_ps = param.hp / (param.hp + param.hs);
% Harmonic mean for the coefficients at the interface with the separator
De_eff_inter_ps = (De_eff_p(end)*De_eff_s(1)) / (betaD_ps*De_eff_s(1)+(1-betaD_ps)*De_eff_p(end));
% De_eff_p at CV edges in positive electrode
De_eff_p_inter = [De_eff_p_inter;De_eff_inter_ps];

%%  Diffusion coefficients for the separator

% Beta coefficient within the separator
betaD_s = 0.5;
% Harmonic mean for the coefficients within the separator
De_eff_s_inter = (De_eff_s(1:end-1).*De_eff_s(2:end)) ./ (betaD_s*De_eff_s(2:end)+(1-betaD_s)*De_eff_s(1:end-1));
% Beta coefficient for the interface with the negative electrode
betaD_sn = param.hs / (param.hs + param.hn);
% Harmonic mean for the coefficients at the interface with the negative electrode
De_eff_inter_sn = (De_eff_s(end)*De_eff_n(1)) / (betaD_sn*De_eff_n(1)+(1-betaD_sn)*De_eff_s(end));
% De_eff_s at CV edges in separator
De_eff_s_inter = [De_eff_s_inter;De_eff_inter_sn];

%% Diffusion coefficients for the negative electrode

% Beta coefficient within the negative electrode
betaD_n  = 0.5;
% De_eff_p at CV edges in negative electrode
De_eff_n_inter = (De_eff_n(1:end-1).*De_eff_n(2:end)) ./ (betaD_n*De_eff_n(2:end)+(1-betaD_n)*De_eff_n(1:end-1));

end

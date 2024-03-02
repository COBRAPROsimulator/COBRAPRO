function t1 = t1_Nyman(ce, param)
% t1_Nyman evaluates the transference number as a function of electrolyte
% concentration (ce). Function from [1].
%
% Reference:
% [1] A. Nyman, M. Behm, and G. Lindbergh, "Electrochemical characterisation and modelling of the mass transport phenomena in LiPF6–EC–EMC electrolyte,” Electrochimica Acta, vol. 53, no. 22, pp. 6356–6365, Sep. 2008, doi: 10.1016/j.electacta.2008.04.023.
%
%   Inputs:
%       ce: normalized electrolyte concentration [-]
%       param: parameter structure
%
%   Output:
%       t1: transference number [-]

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

% Convert to [mol.m-3]
ce = ce.*param.c0;

t1 = -0.1287*(ce/1000).^3+0.4106*(ce/1000).^2-0.4717*(ce/1000)+0.4492;

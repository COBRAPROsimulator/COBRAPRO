function t1 = t1_Kremer(ce, param)
% t1_Kremer evaluates the transference number as a function of electrolyte
% concentration (ce). Function from [1].
%
% Reference:
% [1] L. S. Kremer et al., "Manufacturing Process for Improved Ultra‐Thick Cathodes in High‐Energy Lithium‐Ion Batteries,” Energy Technol., vol. 8, no. 2, p. 1900167, Feb. 2020, doi: 10.1002/ente.201900167.
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

numerator = 3.2093e-4*(ce/1e6).^2-1.2798e-6*(ce/1e6)+1.374e-9;
denominator = (ce/1e6).^3-0.0026*(ce/1e6).^2+6.7273e-8*(ce/1e6)+3.17779e-9;

t1 = numerator./denominator;

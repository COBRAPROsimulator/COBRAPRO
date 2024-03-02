function Deff = De_eff_Kremer(ce, param, domain)
% De_eff_Kremer evaluates the effective electrolyte diffusivity as a
% function of electrolyte concentration (ce). Function from [1].
%
% Reference:
% [1] L. S. Kremer et al., "Manufacturing Process for Improved Ultra‐Thick Cathodes in High‐Energy Lithium‐Ion Batteries,” Energy Technol., vol. 8, no. 2, p. 1900167, Feb. 2020, doi: 10.1002/ente.201900167.
%
%   Inputs:
%       ce: normalized electrolyte concentration [-]
%       param: parameter structure
%       domain: string letter indicating battery domain
%           -> 'p': positive electrode
%           -> 's': separator
%           -> 'n': negative electrode
%
%   Output:
%       Deff: effective electrolyte diffusivity [m2.s-1]

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

switch(domain)
    case'p'
    Deff = param.ep^param.brugp *(-1.2564e5*(ce/1e6).^4+1.179e3*(ce/1e6).^3-3.5799*(ce/1e6).^2+0.0027*(ce/1e6)+3.5884e-6)*0.0001;
    case 's'
    Deff = param.es^param.brugs *(-1.2564e5*(ce/1e6).^4+1.179e3*(ce/1e6).^3-3.5799*(ce/1e6).^2+0.0027*(ce/1e6)+3.5884e-6)*0.0001;
    case 'n'
    Deff = param.en^param.brugn *(-1.2564e5*(ce/1e6).^4+1.179e3*(ce/1e6).^3-3.5799*(ce/1e6).^2+0.0027*(ce/1e6)+3.5884e-6)*0.0001;
end


function Kappa_eff = Kappa_eff_Kremer(ce,param,domain)
% Kappa_eff_Kremer evaluates the effective electrolyte conductivity as a
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
%       Kappa_eff: effective electrolyte conductivity [S.m-1]

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
    Kappa_eff = param.ep^param.brugp *(1.0183e11*(ce/1e6).^5-1.4609e9*(ce/1e6).^4+8.3577e6*(ce/1e6).^3-2.2877e4*(ce/1e6).^2+25.3300*(ce/1e6)+5.2392e-5)*100;
    case 's'
    Kappa_eff = param.es^param.brugs *(1.0183e11*(ce/1e6).^5-1.4609e9*(ce/1e6).^4+8.3577e6*(ce/1e6).^3-2.2877e4*(ce/1e6).^2+25.3300*(ce/1e6)+5.2392e-5)*100;
    case 'n'
    Kappa_eff = param.en^param.brugn *(1.0183e11*(ce/1e6).^5-1.4609e9*(ce/1e6).^4+8.3577e6*(ce/1e6).^3-2.2877e4*(ce/1e6).^2+25.3300*(ce/1e6)+5.2392e-5)*100;
end

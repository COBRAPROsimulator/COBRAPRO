function Kappa_eff = Kappa_eff_Chen2020(ce,param,domain)
% Kappa_eff_Chen2020 evaluates the effective electrolyte conductivity as a
% function of electrolyte concentration (ce). Function from [1].
%
% Reference:
% [1] C.-H. Chen, F. Brosa Planella, K. O’Regan, D. Gastol, W. D. Widanage, and E. Kendrick, “Development of Experimental Techniques for Parameterization of Multi-scale Lithium-ion Battery Models,” J. Electrochem. Soc., vol. 167, no. 8, p. 080534, Jan. 2020, doi: 10.1149/1945-7111/ab9050.
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
    Kappa_eff = param.ep^param.brugp *(0.1297*(ce/1000).^3-2.51*(ce/1000).^1.5+3.329*(ce/1000));
    case 's'
    Kappa_eff = param.es^param.brugs *(0.1297*(ce/1000).^3-2.51*(ce/1000).^1.5+3.329*(ce/1000));
    case 'n'
    Kappa_eff = param.en^param.brugn *(0.1297*(ce/1000).^3-2.51*(ce/1000).^1.5+3.329*(ce/1000));
end


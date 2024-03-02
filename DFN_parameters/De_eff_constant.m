function Deff = De_eff_constant(ce, param, domain)
% De_eff_constant evaluates the constant effective electrolyte diffusivity.
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

switch(domain)
    case 'p'
        Deff = param.ep^param.brugp.*param.De*ce./ce;
    case 's'
        Deff = param.es^param.brugs.*param.De*ce./ce;
    case 'n'
        Deff = param.en^param.brugn.*param.De*ce./ce;
end

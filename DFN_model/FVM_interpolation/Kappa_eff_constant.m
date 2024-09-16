function Kappa_eff = Kappa_eff_constant(ce,param,domain)
% Kappa_eff_constant evaluates the constant effective electrolyte conductivity
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

switch(domain)
    case'p'
    Kappa_eff=param.ep^param.brugp.*param.Kappa*ce./ce;
    case 's'
    Kappa_eff=param.es^param.brugs.*param.Kappa*ce./ce;
    case 'n'
    Kappa_eff=param.en^param.brugn.*param.Kappa*ce./ce;
end


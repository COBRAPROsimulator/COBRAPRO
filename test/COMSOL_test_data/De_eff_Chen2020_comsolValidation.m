function Deff = De_eff_Chen2020_comsolValidation(ce, param, domain)
% De_eff_Chen2020_comsolValidation evaluates the effective electrolyte diffusivity as a
% function of electrolyte concentration (ce). Function from [1].
% Used in test_2_comsolValidation.m script to check the code against COMSOL
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
%       Deff: effective electrolyte diffusivity [m2.s-1]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COBRAPRO: Co-simulation Battery Modeling for Accelerated Parameter Optimization

% Copyright (c) 2024 CO-simulation BatteRy modeling for Accelerated PaRameter Optimization (COBRAPRO)
% COBRAPRO is freely distributed software under the MIT License 
% v1.0.0.: Released 02/2024 

% Written by Sara Ha (sungyeon.sara.ha@stanford.edu)
% PI: Prof. Simona Onori (sonori@stanford.edu)
% Stanford Energy Control Group, Stanford University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert to [mol.m-3]
ce = ce.*param.c0;

switch(domain)
    case 'p'
        Deff = param.ep^param.brugp.*(8.794e-11.*(ce/1000).^2-3.972e-10.*(ce/1000)+4.862e-10);
    case 's'
        Deff = param.es^param.brugs.*(8.794e-11.*(ce/1000).^2-3.972e-10.*(ce/1000)+4.862e-10);
    case 'n'
        Deff = param.en^param.brugn.*(8.794e-11.*(ce/1000).^2-3.972e-10.*(ce/1000)+4.862e-10);
end


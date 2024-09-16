function Deff = De_eff_Nyman2008(ce, param, domain)
% De_eff_Nyman2008 evaluates the effective electrolyte diffusivity as a
% function of electrolyte concentration (ce). Function from [1] and used in [2].
%
% Reference:
% [1] A. Nyman, M. Behm, and G. Lindbergh, "Electrochemical characterisation and modelling of the mass transport phenomena in LiPF6–EC–EMC electrolyte,” Electrochimica Acta, vol. 53, no. 22, pp. 6356–6365, Sep. 2008, doi: 10.1016/j.electacta.2008.04.023.
% [2] C.-H. Chen, F. Brosa Planella, K. O’Regan, D. Gastol, W. D. Widanage, and E. Kendrick, “Development of Experimental Techniques for Parameterization of Multi-scale Lithium-ion Battery Models,” J. Electrochem. Soc., vol. 167, no. 8, p. 080534, Jan. 2020, doi: 10.1149/1945-7111/ab9050.
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

% Function valid for ce = 500 to 2000 [mol.m-3] 
if isnumeric(ce)
    if any(ce < 500)
        ce(ce < 500) = 500;
    elseif any(ce > 2000)
        ce(ce > 2000) = 2000;
    end
end

switch(domain)
    case 'p'
        Deff = param.ep^param.brugp.*(8.794e-11.*(ce/1000).^2-3.972e-10.*(ce/1000)+4.862e-10);
    case 's'
        Deff = param.es^param.brugs.*(8.794e-11.*(ce/1000).^2-3.972e-10.*(ce/1000)+4.862e-10);
    case 'n'
        Deff = param.en^param.brugn.*(8.794e-11.*(ce/1000).^2-3.972e-10.*(ce/1000)+4.862e-10);
end


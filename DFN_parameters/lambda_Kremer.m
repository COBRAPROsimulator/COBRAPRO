function lambda = lambda_Kremer(ce, param)
% lambda_Kremer evaluates the activity coefficient as a function of electrolyte
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
%       lambda: activity coefficient [-]

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

lambda = 4.4825e10*(ce/1e6).^4-4.6529e8*(ce/1e6).^3+1.3801e6*(ce/1e6).^2+89.0407*(ce/1e6)+0.9378;

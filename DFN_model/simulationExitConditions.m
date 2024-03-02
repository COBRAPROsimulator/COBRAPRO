function [exit_indicator, exit_value, exit_message] = simulationExitConditions(voltage, csp_surf, csn_surf, param, curr_dens)
%   simulationExitConditions returns the simulation exit indicator and exit value.
%
%   Inputs:
%       voltage: Voltage vector [V]
%       csp: Normalized concentration in positive electrode particles [-]
%       csn: Normalized concentration in positive electrode particles [-]
%       csp_surf: Normalized surface concentration in positive electrode particles [-]
%       csn_surf: Normalized surface concentration in positive electrode particles [-]
%       param: Parameter structure
%       curr_dens: Current density vector [A.m-2]
%
%   Output:
%       exit_indicator: String indiciating reason for simulation ending
%           -> " ": no exit conditions triggered
%           -> "V_cutoffUpper": Upper cutoff voltage triggered
%           -> "V_cutoffLower": Lower cutoff voltage triggered
%           -> "cspSurf_aboveCspMax: csp_surf greater than csp_max (results in imaginary exchange current density, jp0)
%           -> "csnSurf_aboveCsnMax: csn_surf greater than csn_max (results in imaginary exchange current density, jn0)
%       exit_value: Value corresponding to when the exit condition is triggered
%       exit_message: Message explaining type of error

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

exit_indicator = {};
exit_value = [];
exit_message = {};

% Check if cs_surf > cs_max (square root in BV equation will be negative)
if (csp_surf(end,:) > param.ctp)
    exit_message = [exit_message;{'Positive surface particle concentration above csp_max.\n'}];
    exit_indicator = [exit_indicator;{'cspSurf_aboveCspMax'}];
    exit_value = [exit_value;csp_surf(end,:)];
elseif (csn_surf(end,:) > param.ctn) 
    exit_message = [exit_message;{'Negative surface particle concentration above csn_max.\n'}];
    exit_indicator = [exit_indicator;{'csnSurf_aboveCsnMax'}];
    exit_value = [exit_value;csn_surf(end,:)];
end 

if curr_dens(end) >= 0 % In charging
    if(voltage(end) > param.upperCutoffVoltage)
        exit_message = [exit_message;{'Cell above its upper cut-off voltage.\n'}];
        exit_indicator = [exit_indicator;{'V_cutoffUpper'}];
        exit_value = [exit_value;voltage(end)];
    end
elseif curr_dens(end) <= 0 % In discharging
    if(voltage(end) < param.lowerCutoffVoltage)
        exit_message = [exit_message;{'Cell below its lower cut-off voltage.\n'}];
        exit_indicator = [exit_indicator;{'V_cutoffLower'}];
        exit_value = [exit_value;voltage(end)];
    end
end
 
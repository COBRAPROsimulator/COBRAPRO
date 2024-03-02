% -------------------------------------------------------------------------
% test_1_casadiCheck.m
% -------------------------------------------------------------------------
%   test_1_casadiCheck determines if CasADi is working properly. This script 
%   automatically checks if CasADi is installed correctly.

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

import casadi.*

try
    x = SX.sym('x',1,3);
    fprintf('\ntest_1 successful: CasADi is working properly!\n\n');
catch error
    fprintf('\ntest_2 unsuccesful.\n')
    fprintf(['\nThere was an error:\n' error.message '\n\n']);
    fprintf(['The error identifier is:\n' error.identifier '\n\n']);
end

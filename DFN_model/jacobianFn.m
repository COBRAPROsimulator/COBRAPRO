function [J, flag, new_data] = jacobianFn(t, x, xp, ~, cj, data)
%	jacobianFn calculates the Jacobian matrix of the DFN model
%
%   Inputs:
%       t: time 
%       x: Vector consisting of algebraic and differential variables
%       xp: Vector consisting of derivative of algebraic and differential variables
%       cj: required by IDA to compute Jacobian matrix
%       data: required by IDA to extract Jocobian function
%
%   Output:
%       J: Jacobian matrix
%       flag: Required by IDA solver but not used in code
%       new_data: Required by IDA solver but not used in code

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

% Obtaini Jacobian function from data
fJ = data.fJ;

% Compute the Jacobian for the given x and xp at time t
J  = full(fJ(t,x,xp,cj));

% Required by IDA solver but not used in this code
flag        = 0;
new_data    = [];

end
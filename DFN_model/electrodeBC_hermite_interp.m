function states_cc = electrodeBC_hermite_interp(states,curr_dens,param,domain)
%   electrodeBC_hermite_interp computes the value at the
%   CC|positive or negative|CC using 3rd order hermite interpolation.
%
%   Inputs:
%       states: vector of states in positive or negative electrode 
%       curr_dens: current density [A.m-2]
%       param: parameter structure
%       domain: string letter indicating electrode domain
%           -> 'p': positive electrode
%           -> 'n': negative electrode
%
%   Output:
%       states_cc: vector of states at CC|positive or negative|CC 

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

switch domain
    case 'p'
        x_1 = param.hp*(0.5); % x distance from i = 0.5 to i = 1
        x_2 = param.hp*(1.5); % x distance from i = 0.5 to i = 2
        p1 = (1+2*(-x_2)/(x_1-x_2))*states(:,2)+(-x_2)*((states(:,3)-states(:,1))/(2*param.hp));
        p2 = ((-x_1)/(x_2-x_1))^2;
        p3 = (1+2*(x_1)/(x_1-x_2))*states(:,1)+(-x_1)*(-curr_dens/param.sigma_effp);
        p4 = ((-x_2)/(x_1-x_2))^2;
        states_cc = p1*p2 + p3*p4;
    case 'n'
        x_Nn = param.hn*(param.Nn-0.5); % x distance from i = 0.5 to i = Nn
        x_Nn_1 = param.hn*(param.Nn-1.5); % x distance from i = 0.5 to i = Nn-1
        n1 = (1+2*(param.ln-x_Nn_1)/(x_Nn-x_Nn_1))*states(:,param.Nn-1)+(param.ln-x_Nn_1)*((states(:,param.Nn)-states(:,param.Nn-2))/(2*param.hn));
        n2 = ((param.ln-x_Nn)/(x_Nn_1-x_Nn))^2;
        n3 = (1+2*(x_Nn-param.ln)/(x_Nn-x_Nn_1))*states(:,param.Nn)+(param.ln-x_Nn)*(-curr_dens/param.sigma_effn);
        n4 = ((param.ln-x_Nn_1)/(x_Nn-x_Nn_1))^2;
        states_cc = n1*n2 + n3*n4; 
end
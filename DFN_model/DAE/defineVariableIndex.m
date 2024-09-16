function param = defineVariableIndex(param)
%   defineVariableIndex defines index and length of differential and algebraic 
%   variables and stores in param structure.
%
%   Inputs:
%       param: Parameter structure
%
%   Output:
%       param: Parameter structure

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

% Extract discretized node number.
Np=param.Np;
Ns=param.Ns;
Nn=param.Nn;
Nrp=param.Nrp;
Nrn=param.Nrn;

% Calculate variable length for each variable type.
ce_length = Np+Ns+Nn;     
csp_length = Np*Nrp;         
csn_length = Nn*Nrn;       
phisp_length = Np;
phisn_length = Nn;
phie_length = Np+Ns+Nn;   
jp_length = Np;
jn_length = Nn;
I_density_length = 1;

% Store variable length in param.
param.ce_length = ce_length;
param.csp_length = csp_length;
param.csn_length = csn_length;
param.phisp_length = phisp_length;
param.phisn_length = phisn_length;
param.phie_length = phie_length;
param.jp_length = jp_length;
param.jn_length = jn_length;
param.jn_length = jn_length;
param.I_density_length = I_density_length;

% Calculate last index of each variable.
ce_last_ind = ce_length;
csp_last_ind = ce_last_ind + csp_length;
csn_last_ind = csp_last_ind + csn_length;
phisp_last_ind = csn_last_ind + phisp_length;
phisn_last_ind = phisp_last_ind + phisn_length;
phie_last_ind = phisn_last_ind + phie_length;
jp_last_ind = phie_last_ind + jp_length;
jn_last_ind = jp_last_ind + jn_length;
I_density_ind = jn_last_ind + I_density_length;

% Store position of differential and algebraic variables in param.
param.ce_ind = 1:ce_last_ind;
param.csp_ind = ce_last_ind+1:csp_last_ind;
param.csn_ind = csp_last_ind+1:csn_last_ind;
param.phisp_ind = csn_last_ind+1:phisp_last_ind;
param.phisn_ind = phisp_last_ind+1:phisn_last_ind;
param.phie_ind = phisn_last_ind+1:phie_last_ind;
param.jp_ind = phie_last_ind+1:jp_last_ind;
param.jn_ind = jp_last_ind+1:jn_last_ind;
param.I_density_ind = I_density_ind;

% Store total number of variables in param.
total_variables = I_density_ind;
param.total_variables = total_variables;

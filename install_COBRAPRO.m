function install_COBRAPRO
%   install_COBRAPRO installs SUNDIALS 2.6.2 and adds all the necessary
%   folders to the MATLAB path. Before running this function, the 
%   SUNDIALS 2.6.2 and CasADi folders should be inside same directory 
%   as install_COBRAPRO.m

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

% SUNDIALS 2.6.2 folder name
sundials_folder_name = 'sundials-2.6.2';

% Install SUNDIALS 2.6.2
%--------------------------------------------------------------------------
% In your Command Window you will see the following prompts:

% MEX files will be compiled and built using the above options
%   Proceed? (y/n)-> type "y"

% Compile CVODES interface? (y/n) -> type "n"
% Compile IDAS interface? (y/n) -> type "y"
% Compile KINSOL interface? (y/n) -> type "n"

% MEX files were successfully created.
%   Install toolbox? (y/n) -> type "y"

% Enter return to cancel the installation.
%   Installation directory: -> just hit enter
%--------------------------------------------------------------------------

if isempty(genpath(sundials_folder_name))
    fprintf('Failed: SUNDIALS 2.6.2 folder not located inside COBRAPRO folder.\n')
    return
else
    run sundials-2.6.2/sundialsTB/install_STB.m
    addpath(genpath(sundials_folder_name)) 
end
fprintf('SUNDIALS installed sucessfully!\n')

% Add CasADi folder to path
fprintf('\n\n Setting up CasADi...\n')
casadi_folder_name = input('\nType the name of your CasADi folder (case-sensitive): ','s');
if isempty(genpath(casadi_folder_name))
    fprintf('Failed: CasADi folder not located inside COBRAPRO folder.\n')
    return
else
    addpath(genpath(casadi_folder_name))
end

% All necessary folders to path
folder = fileparts(which('install_COBRAPRO.m'));
addpath(genpath(folder));
% Save path
savepath

fprintf('\n\nCOBRAPRO installed successfully! To get started, try running scripts in the Examples folder.\n\n')

end

function param_LSA_unCorr = unCorr_parameters(param_LSA, S_V_SOC, corr_V_SOC_matrix, beta)
%   unCorr_parameters starts with the parameter with the highest sensitivty and
%   removes the parameters in param_LSA if that parameter has a correlation 
%   coefficient higher than the threshold (beta). This process repeats
%   until all the parameters have been investigated in param_LSA.
%
%   Inputs:
%       param_LSA: cell containing names of the parameters
%       S_V_SOC: vector consisting of sensitivity values for each parameter
%       corr_V_SOC_matrix: absolute value of the correlation matrix
%       beta: correlation coefficient threshold [-]
%
%   Output:
%       param_LSA_unCorr: cell containing names of the uncorrelated parameters

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

    % Initialize param_LSA_unCorr with param_LSA
    param_LSA_unCorr = param_LSA;

    % Sort param_LSA_unCorr based on S_V_SOC from highest to lowest sensitivity
    [~, sortedIndices] = sort(S_V_SOC, 'descend');
    param_LSA_unCorr_sorted = param_LSA_unCorr(sortedIndices);

    % Define a cell that will keep track of the remaining parameters when the 
    % parameters are eliminated from the cell
    theta_remaining = param_LSA_unCorr_sorted;

    % Initialize a logical array to keep track of the indices to remove for theta_remaining (length changing during while loop) 
    toKeep_remove = true(size(param_LSA_unCorr_sorted));
    % Initialize a logical array to keep track of the indices to remove for param_LSA_unCorr
    toKeep = true(size(param_LSA_unCorr_sorted));

    % Determine uncorrelated parameters
    while length(theta_remaining) > 1
        % Find actual index for the first parameter being investigated (will always be the first index of theta_remaining by construction) 
        i = strcmp(param_LSA_unCorr,theta_remaining{1});
        % Compare with each parameter in the list of remaining parameters
        for j = 2:length(theta_remaining)
            % Find actual index for the second parameter being investigated
            k = strcmp(param_LSA_unCorr,theta_remaining{j});
            % Extract the correlation coefficient for the given two parameters
            C_value = corr_V_SOC_matrix(i,k);
            % If correlation is above threshold, mark it for removal
            if C_value > beta
                toKeep(k) = false;
                toKeep_remove(j) = false;
            end
        end
    % Remove the first element in the toKeep_remove vector since we are done investigating the first parameter in this vector 
    toKeep_remove(1) = false;
    % Update theta_remaining to remove the parameter that was investigated
    % and the parameters correlated with that parameter
    theta_remaining = theta_remaining(toKeep_remove);
    % Re-initialize logical array to match new size of theta_remaining
    toKeep_remove = true(size(theta_remaining));
    end

    % Output uncorrelated parameters
    param_LSA_unCorr = param_LSA_unCorr(toKeep);
end

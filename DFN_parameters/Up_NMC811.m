function Up = Up_NMC811(theta)
% Up_NMC811 is the experimentally determined open-circuit potential (OCP)
% of the NMC811 cathode from [1]. OCP was determined from half-cell GITT
% tests in [1].
%
% Reference:
% [1] C.-H. Chen, F. Brosa Planella, K. O’Regan, D. Gastol, W. D. Widanage, and E. Kendrick, “Development of Experimental Techniques for Parameterization of Multi-scale Lithium-ion Battery Models,” J. Electrochem. Soc., vol. 167, no. 8, p. 080534, Jan. 2020, doi: 10.1149/1945-7111/ab9050.
%
%   Inputs:
%       theta: normalized surface particle concentration in cathode [-]
%
%   Output:
%       Up: open-circuit potential in cathode [V]

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

Up=-0.8090.*theta+4.4875-0.0428.*tanh(18.5138.*(theta-0.5542))-17.7326.*tanh(15.7890.*(theta-0.3117))+17.5842.*tanh(15.9308.*(theta-0.3120));

end
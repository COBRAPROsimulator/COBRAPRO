function Un = Un_Graphite_Si_Chen2020(theta)
% Un_Graphite_Si is the experimentally determined open-circuit potential (OCP)
% of the Graphite-Silicon anode from [1]. OCP was determined from half-cell GITT
% tests in [1].
%
% Reference:
% [1] C.-H. Chen, F. Brosa Planella, K. O’Regan, D. Gastol, W. D. Widanage, and E. Kendrick, “Development of Experimental Techniques for Parameterization of Multi-scale Lithium-ion Battery Models,” J. Electrochem. Soc., vol. 167, no. 8, p. 080534, Jan. 2020, doi: 10.1149/1945-7111/ab9050.
%
%   Inputs:
%       theta: normalized surface particle concentration in anode [-]
%
%   Output:
%       Un: open-circuit potential in anode [V]

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

Un=1.9793.*exp(-39.3631.*theta)+0.2482-0.0909.*tanh(29.8538.*(theta-0.1234))-0.04478.*tanh(14.9159.*(theta-0.2769))-0.0205.*tanh(30.4444.*(theta-0.6103));

end

function param = Parameters_LG_INR21700_M50
% Parameters_LG_INR21700_M50 loads all the nominal DFN parameters obtained 
% from [1]. The cell is from LG INR21700 M50 its chemistry is 
% NMC811|Graphite-Si with a nominal capacity of 4.85 Ah. This function
% also loads all the DFN simulation settings into the param structure.
%
% References:
% [1] C.-H. Chen, F. Brosa Planella, K. O’Regan, D. Gastol, W. D. Widanage, and E. Kendrick, “Development of Experimental Techniques for Parameterization of Multi-scale Lithium-ion Battery Models,” J. Electrochem. Soc., vol. 167, no. 8, p. 080534, Jan. 2020, doi: 10.1149/1945-7111/ab9050.
% [2] S. Ha, G. Pozzato, and S. Onori, "Electrochemical characterization tools for lithium-ion batteries," J Solid State Electrochem, Nov. 2023, doi: 10.1007/s10008-023-05717-1.
% [3] G. Pozzato, A. Allam, and S. Onori, "Lithium-ion battery aging dataset based on electric vehicle real-driving profiles," Data in Brief, vol. 41, p. 107995, Apr. 2022, doi: 10.1016/j.dib.2022.107995.
%
%   Inputs:
%
%   Output:
%       param: structure containing nominal parameters of LG 21700 M50
%       cells and the DFN simulation settings

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

%--------------------------------------------------------------------------
% General Parameters
%--------------------------------------------------------------------------
% Faraday constant [C/mol]
param.F=96487; 
% Ideal gas constant [J/(mol.K)]
param.R=8.3143;    
% Temperature [K]
param.T=296.15; % 23 Celsius

%--------------------------------------------------------------------------
% Open-circuit voltage function: U(cs_surf/cs_max) [V] [Chen 2020]
%--------------------------------------------------------------------------
param.Up = @Up_NMC811_Chen2020; 
param.Un = @Un_Graphite_Si_Chen2020; 

%--------------------------------------------------------------------------
% Cut-off voltage [V] [Chen 2020]
%--------------------------------------------------------------------------
param.lowerCutoffVoltage = 2.5; 
param.upperCutoffVoltage = 4.2; 

%--------------------------------------------------------------------------
% Stoichiometric parameters [-] [Chen 2020]
%--------------------------------------------------------------------------
% Stoichiometric limit in cathode at 100% SOC
param.theta100_p = 0.27; 
% Stoichiometric limit in anode at 100% SOC
param.theta100_n = 0.9014; 
% Stoichiometric limit in cathode at 0% SOC
param.theta0_p = 0.9084; 
% Stoichiometric limit in anode at 0% SOC
param.theta0_n = 0.0279; 

%--------------------------------------------------------------------------
% Contact resistance [Ohms] [Ha 2023, Pozzato 2022]
%--------------------------------------------------------------------------
% Obtained from EIS characterization from cell W8 (LG INR21700 M50T) at Diag.#1 (fresh cell) [Ha 2023, Pozzato 2022]
param.Rc = 0.02339; 

%--------------------------------------------------------------------------
% Design Parameters [Chen 2020]
%--------------------------------------------------------------------------
% Filler fraction at cathode [-]
param.ep_f = 0;      
% Filler fraction at anode [-]
param.en_f = 0;              
% Porosity at cathode [-]
param.ep = 0.335; 
% Porosity at separator [-]
param.es = 0.47;        
% Porosity at anode [-]
param.en = 0.25;            
% Cathode solid phase volume fraction [-] 
param.ep_s = 1-param.ep-param.ep_f;  
% Anode solid phase volume fraction [-] 
param.en_s = 1-param.en;                
% Bruggeman coefficient at cathode [-]
param.brugp = 1.5;   
% Bruggeman coefficient at separator [-]
param.brugs = 1.5;   
% Bruggeman coefficient at anode [-]
param.brugn = 1.5;           
% Thickness at cathode [m]
param.lp = 75.6e-6;   
% Thickness at separator [m]
param.ls = 12e-6;       
% Thickness at anode [m]
param.ln = 85.2e-6;             
% Radius of solid particle at cathode [m] 
param.Rpp = 5.22e-6;               
% Radius of solid particle at anode [m]
param.Rpn = 5.86e-6;               
% Particle surface area to volume at cathode [m-1]
param.ap = (3/param.Rpp)*param.ep_s;      
% Particle surface area to volume at anode [m-1]
param.an = (3/param.Rpn)*param.en_s;    
% Cell cross-sectional area [m2]
param.Acell = (0.1027 + 0.10465)/2;           

%--------------------------------------------------------------------------
% Cell nominal capacity [Ah]
%--------------------------------------------------------------------------
% Nominal capacity can be determined from:
% Option 1. Calculating discharge capacity (Q_dis) from C/20 capacity test of a fresh cell [Ha 2023]
% Option 2. If Q_dis not available, you may also use the nominal capacity reported by the cell manufacturer
%--------------------------------------------------------------------------
% Calculated from C/20 capacity test for cell W8 at Diag.#1 [Ha 2023]    
param.Q_nom = 4.8768;     

%--------------------------------------------------------------------------
% Concentration parameters [Chen 2020] 
%--------------------------------------------------------------------------    
% Initial electrolyte concentration [mol.m-3]
param.c0 = 1000;       
% Maximum solid phase concentration at cathode [mol.m-3]
param.ctp = 63104;
% Maximum solid phase concentration at anode [mol.m-3]
param.ctn = 33133;

%--------------------------------------------------------------------------
% Transport parameters [Chen 2020] 
%--------------------------------------------------------------------------    
% Solid phase conductivity at cathode [S/m]
param.sigmap = 0.2707;       
% Solid phase conductivity at anode [S/m]
param.sigman = 286.7;       
% Solid phase conductivity at cathode [S/m]
param.sigma_effp = param.ep_s*param.sigmap;
% Solid phase conductivity at anode [S/m]   
param.sigma_effn = param.en_s*param.sigman;

% Solid particle diffusivity at cathode [m2/s]
param.Dsp = 4e-15;
% Solid particle diffusivity at anode [m2/s]
param.Dsn = 3.3e-14;

%--------------------------------------------------------------------------
% Option 1: Functional electrolyte diffusion coefficient: De_eff(ce) [m2/s]
%--------------------------------------------------------------------------
% param.De_eff = @De_eff_Chen2020; 
%--------------------------------------------------------------------------
% Option 2: Constant electrolyte diffusion coefficient: De_eff [m2/s]
%--------------------------------------------------------------------------
param.De = 3.7621e-10;
param.De_eff = @De_eff_constant;

%--------------------------------------------------------------------------
% Option 1: Functional electrolyte conductivity coefficient: kappa_eff(ce) [S/m] [Chen 2020]
%--------------------------------------------------------------------------
% param.Kappa_eff = @Kappa_eff_Chen2020;
%--------------------------------------------------------------------------
% Option 2: Constant electrolyte conductivity coefficient: kappa_eff [S/m] [Chen 2020]
%--------------------------------------------------------------------------
param.Kappa = 0.95;
param.Kappa_eff = @Kappa_eff_constant;

%--------------------------------------------------------------------------
% Option 1: Functional transference number function: t1(ce) [-] [Chen 2020]
%--------------------------------------------------------------------------
% param.t1 = @t1_Kremer;
%--------------------------------------------------------------------------
% Option 2: Constant functional transference number: t1 [-] [Chen 2020]
%--------------------------------------------------------------------------
param.t1_constant = 0.2523;
param.t1= @t1_constant;

%--------------------------------------------------------------------------
% Option 1: Functional activity term: lambda(ce) [-] [Chen 2020]
%--------------------------------------------------------------------------
% param.lambda = @lambda_Kremer; 
%--------------------------------------------------------------------------
% Option 2: Scalar activity term: lambda [-] [Chen 2020]
%--------------------------------------------------------------------------
param.lambda_constant = 1;
param.lambda = @lambda_constant;

%--------------------------------------------------------------------------
% Kinetic parameters [m2.5/(mol0.5.s)] [Chen 2020]
%--------------------------------------------------------------------------
% Reaction rate constant at cathode
param.kp = 3.42e-6/param.F;      
% Reaction rate constant at anode
param.kn = 6.48e-7/param.F;          

%--------------------------------------------------------------------------
% Node numbers
%--------------------------------------------------------------------------
param.Np = 10;
param.Ns = 10;
param.Nn = 10;
param.Nrp = 20;
param.Nrn = 20;

%--------------------------------------------------------------------------
% Enter number of desired workers in a parallel pool
%--------------------------------------------------------------------------
param.pso_workers = 24;

%--------------------------------------------------------------------------
% A single-step iteration-free initialization approach (Patented by Venkat Subramanian Group at UT Austin)- Allows for Only Academic purpose
%--------------------------------------------------------------------------
% Perturbation time [s] (time allowed for perturbation approach to find consistent IC)
param.initime = 3; 
% Perturbation coefficient
param.mu = 10^(-3); 
% Switch coefficient
param.q = 1000; 

%--------------------------------------------------------------------------
% Initialization method
%--------------------------------------------------------------------------
% param.init_Method = 'SS'          -> Single-step method
% param.init_Method = 'IDACalcIC'   -> SUNDIALS IDACalcIC method
%--------------------------------------------------------------------------
param.init_Method = 'SS';

%--------------------------------------------------------------------------
% IDA solver type
%--------------------------------------------------------------------------
% param.IDA_type = 1 -> OneStep (IDA solver determines efficient time step--usually faster than Normal type)
% param.IDA_type = 2 -> Normal  (user defined time step)
%--------------------------------------------------------------------------
param.IDA_type = 1;

%--------------------------------------------------------------------------
% Choose your sampling time for simulation
%--------------------------------------------------------------------------
% -> Used to interpolate the final simulation. This particularly
% useful when using OneStep IDA solver method ("param.IDA_type = 1") since
% the time step is variable when using the OneStep method. Sometimes a large 
% time step during the simulation makes the voltage curve look not smooth. 
% Make the voltage curve "smooth" by using a uniform smaller sampling time 
% to interpolate the final results.
%--------------------------------------------------------------------------
param.deltaT_sim = 1; % [s] 

%--------------------------------------------------------------------------
% Original solution memory 
%--------------------------------------------------------------------------
% param.save_memory = 0  -> Discard original solution (recommend for reduced memory of output solution)
% param.save_memory = 1  -> Save original solution (for single-step method, this will keep the solution during the initialization time)
%--------------------------------------------------------------------------
param.save_memory = 0;

%--------------------------------------------------------------------------
% Solid particle concentration discretization method (csp and csn)
%--------------------------------------------------------------------------
% param.cs_discret = 'FVM' -> finite volume method
% param.cs_discret = 'FDM' -> finite difference method
%--------------------------------------------------------------------------
param.cs_discret = 'FVM';

%--------------------------------------------------------------------------
% Voltage interpolation mode (only for cs_discret = 'FVM')
%--------------------------------------------------------------------------
% param.VinterpMode = 1 -> no interpolation
% param.VinterpMode = 2 -> linear interplation
% param.VinterpMode = 3 -> 3rd order Hermite interpolation (keep at this mode for highest accuracy)
%--------------------------------------------------------------------------
param.VinterpMode = 3;

%--------------------------------------------------------------------------
% Current mode 
%--------------------------------------------------------------------------
% param.CurrentMode = 1     -> constant current
% param.CurrentMode = 2     -> variable current, e.g., UDDS
%--------------------------------------------------------------------------
param.CurrentMode = 1;

%--------------------------------------------------------------------------
% Print simulation time (for effecient PSO, suppress simulating time printing)
%--------------------------------------------------------------------------
% param.sim_print = 0     -> supress simulation time
% param.sim_print = 1     -> print simulation time
%--------------------------------------------------------------------------
param.sim_print = 0;

%--------------------------------------------------------------------------
% Print cost function value (useful if you want to see the J value calculated for every particle)
%--------------------------------------------------------------------------
% param.J_print = 0     -> supress J value during PSO
% param.J_print = 1     -> print J value during PSO
%--------------------------------------------------------------------------
param.J_print = 1;

%--------------------------------------------------------------------------
% IDA tolernace
%--------------------------------------------------------------------------
param.rtol = 1e-6;
param.atol = 1e-9;

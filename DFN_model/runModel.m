function [output,param] = runModel(x_init,param,t0,tf,SOC_init,currentDensity)
%   runModel is the main function that computes and stores the DFN model
%   simulation results in the output structure
%
%   Inputs:
%       x_init: initial condition guess vector 
%           -> if vector is empty, the code calculates the initial guess based on open-circuit voltage conditions
%           -> if vector is not empty, the initial states are calculated from a previous simulation 
%       param: parameter structure
%       t0: initial time [s]
%       tf: final time [s]
%       SOC_init: initial SOC [-]
%       currentDensity: current density input [A.m-2]
%
%   Output:
%       output: Model output structure
%           - output.t: time vector [s] [Mx1]
%           - output.voltage: voltage vector [V] [Mx1]
%           - output.curr_dens: current density vector [A.m-2] [Mx1]
%           - output.soc_p_bulk: state-of-charge in positive electrode [-] [Mx1]
%           - output.soc_n_bulk: state-of-charge in negative electrode [-] [Mx1]
%           - output.ce: electrolyte concentration matrix [mol.m-3] [Mx(Np+Ns+Nn)]
%           - output.csp: positive particle solid concentration matrix [mol.m-3] [Mx(NrpxNp)]
%           - output.csn: negative particle solid concentration matrix [mol.m-3] [Mx(NrnxNn)]
%           - output.csp_surf: positive particle surface solid concentration matrix [mol.m-3] [MxNp]
%           - output.csp_avg: positive particle average solid concentration matrix [mol.m-3] [MxNp]
%           - output.csn_surf: negative particle surface solid concentration matrix [mol.m-3] [MxNn]
%           - output.csn_avg: ngeative particle average solid concentration matrix [mol.m-3] [MxNn]
%           - output.phisp: positive solid potential matrix [V] [MxNp]
%           - output.phisn: negative solid potential matrix [V] [MxNn]
%           - output.phie: electrolyte potential matrix [V] [MxNn]
%           - output.jp: positive pore-wall flux [mol.m-3.s-1] [MxNp]
%           - output.jn: negative pore-wall flux [mol.m-3.s-1] [MxNn]
%           - output.Up: open-circuit potential in positive electrode [V] [MxNp]
%           - output.Un: open-circuit potential in negative electrode [V] [MxNn]
%           - output.etap: positive overpotential [V] [MxNp]
%           - output.etan: negative overpotential [V] [MxNp]
%           -> M: number of discrete points in time
%           -> Np: number of discretization points in positive electrode
%           -> Ns: number of discretization points in separator
%           -> Nn: number of discretization points in negative electrode
%           -> Nrp: number of discretization points in positive particle
%           -> Nrn: number of discretization points in negative particle
%       param: parameter structure

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

tStart_DAE_build = tic;

% Extract the node numbers for cleaner code.
Np=param.Np;
Ns=param.Ns;
Nn=param.Nn;
Nrp=param.Nrp;
Nrn=param.Nrn;

% Define position of differential and algebraic variables and store in
% param structure.
param = defineVariableIndex(param);

%% Step size
hp=param.lp/Np;                 % Step size at positive (FVM)
hs=param.ls/Ns;                 % Step size at separator (FVM)
hn=param.ln/Nn;                 % Step size at negative (FVM)
rp=param.Rpp/Nrp;               % Step size at positive particle (FVM)
rn=param.Rpn/Nrn;               % Step size at negative particle (FVM)

% Store in parameter structure
param.hp=hp;
param.hs=hs;
param.hn=hn;
param.rp=rp;
param.rn=rn;

%% Spatial points
param.x_intval_p = hp.*((1:Np) - 0.5);
param.x_intval_s = param.lp + hs.*((1:Ns) - 0.5);
param.x_intval_n = param.lp + param.ls + hn.*((1:Nn) - 0.5);
param.x_intval_psn = [param.x_intval_p, param.x_intval_s, param.x_intval_n];
if strcmp(param.cs_discret,'FVM')
    param.r_intval_p = rp.*((1:Nrp) - 0.5);
    param.r_intval_n = rn.*((1:Nrn) - 0.5); 
elseif strcmp(param.cs_discret,'FDM')
    param.r_intval_p = rp.*((1:Nrp+1) - 1);
    param.r_intval_n = rn.*((1:Nrn+1) - 1); 
end

%% Define I_density in param
switch(param.CurrentMode)
    case 1 % constant current
        param.currentDensity = @(t,param)currentDensity;
    case 2 % variable current
        param.currentDensity = param.currentDensityFunction;
        % for single-step approach, need to account for initialization time in the time varying current
        if strcmp(param.init_Method,'SS')
            % If 'SS' make current during param.initime the same as I_data(1)
            param.I_data = [ones(param.initime,1)*param.I_data(1); param.I_data];
            param.t_data = [(0:(param.initime-1))'; param.t_data + param.initime];
        end
end

%% Define symbolic variables used in IDA solver
import casadi.*

t_SX    = SX.sym('t');                           % time
x_SX    = SX.sym('x',param.total_variables);     % variables (differential + algebraic)
xp_SX   = SX.sym('xp',param.total_variables);    % derivative of variables (differential + algebraic)
cj      = SX.sym('cj');                          % required for Jacobian calculation

%% Extract initialization method 
% Store init_Method in ida_user_data for IDA solver
ida_user_data.init_Method = param.init_Method;

%% Get number of algebraic and differential variables for IDACalcIC

if strcmp(param.init_Method,'IDACalcIC')
    % Define algebraic and differential variables 
    %   id:1-> differential variables (according to SUNDIALS TB manual)
    %   id:0-> algebraic variables (according to SUNDIALS TB manual)
    id = zeros(param.total_variables,1);
    id(1:(param.ce_length+param.csp_length+param.csn_length)) = 1;
end

%% Store symbolic implicit ode function in UserData (Only for single-step method)
% Create implicit ode by perturbing algebraic equations!

if strcmp(param.init_Method,'SS')
    % 1. Call function to calculate symbolic algebraic equations
    z_tot = alg_res(x_SX, param);
    
    % 2. Residual of symbolic implicit ode by perturbing algebraic equations
    res_perAlg = param.mu*(jacobian(z_tot,x_SX)*xp_SX)+z_tot;
    
    % 3. Define a casadi Function to evaluate implicit ode given set of algebraic variables
    f_res_perAlg = Function('f_res_perAlg',{x_SX,xp_SX},{res_perAlg});
    
    % 4. Store function into user data so that IDA will use it for evaluation of the overall residual in alg_res.m
    ida_user_data.res_perAlg = f_res_perAlg;
end

%% Create Jacobian matrix in symbolic form and store in UserData
% Store param in ida_user_data for IDA solver
ida_user_data.param = param;

% 1. Calculate the symbolic DAE residual equations
% dx_tot = F(t,x,xp) where F(t,x,xp)=0 is the DAE written in implicit form
[dx_tot, ~, ~] = diff_res(t_SX,x_SX,xp_SX,ida_user_data);

% 2. Compute the Jacobian matrix
J = jacobian(dx_tot,x_SX) + cj*jacobian(dx_tot,xp_SX);

% 3. Define "JacFun" function, which computes the Jacobian matrix
JacFun = Function('fJ',{t_SX,x_SX,xp_SX,cj},{J});

% 4. Store JacFun in ida_user_data. IDA solver will call JacFun to compute the Jacobian matrix
ida_user_data.fJ = JacFun;

%% Initial condition guess

% IC guess is computed using the OCV conditions
if isempty(x_init) 
    % Electrolyte concentration  
    x0(param.ce_ind,1)=1;             
    % Positive concentration (based on SOC_init)
    theta_p_init = SOC_init*(param.theta100_p-param.theta0_p)+param.theta0_p;
    x0(param.csp_ind,1)= theta_p_init;        
    % Negative concentration (based on SOC_init)       
    theta_n_init = SOC_init*(param.theta100_n-param.theta0_n)+param.theta0_n;
    x0(param.csn_ind,1)=theta_n_init; 
    % Positive solid phase potential 
    x0(param.phisp_ind,1)= param.Up(theta_p_init)-param.Un(theta_n_init);        
    % Negative solid phase potential
    x0(param.phisn_ind,1)=0;     
    % Liquid potential
    x0(param.phie_ind,1)=0;  
    % Pore-wall flux
    x0(param.jp_ind,1)=0; 
    x0(param.jn_ind,1)=0;
    % Current density
    x0(param.I_density_ind,1)= param.currentDensity(0,param);
else % Initial conditions are carried over from a previous simulation 
    % Electrolyte concentration  
    x0(param.ce_ind,1)=x_init(param.ce_ind);             
    % Positive concentration 
    x0(param.csp_ind,1)= x_init(param.csp_ind);        
    % Negative concentration      
    x0(param.csn_ind,1)=x_init(param.csn_ind); 
    % Positive solid phase potential 
    x0(param.phisp_ind,1)= x_init(param.phisp_ind);        
    % Negative solid phase potential
    x0(param.phisn_ind,1)=x_init(param.phisn_ind);    
    % Liquid potential
    x0(param.phie_ind,1)=x_init(param.phie_ind);
    % Pore-wall flux
    x0(param.jp_ind,1)=x_init(param.jp_ind);
    x0(param.jn_ind,1)=x_init(param.jn_ind);    
    % Current density
    % x0(param.I_density_ind,1)= x_init(param.I_density_ind);
    x0(param.I_density_ind,1)= param.currentDensity(0,param);
end

% Initialize the xp variables (dx/dt(t=0)=0). 
xp0 = zeros(param.total_variables,1);

%% Setup IDA Solver

if strcmp(param.init_Method,'SS')
    options_IDA = IDASetOptions('RelTol', param.rtol,...
                            'AbsTol', param.atol,...
                            'MaxNumSteps', 1500,...
                            'InitialStep', 1e-3,...
                            'UserData', ida_user_data,...
                            'JacobianFn',@jacobianFn);
elseif strcmp(param.init_Method,'IDACalcIC')
    options_IDA = IDASetOptions('RelTol', param.rtol,...
                            'AbsTol', param.atol,...
                            'MaxNumSteps', 1500,...
                            'InitialStep', 1e-3,...
                            'VariableTypes', id,...
                            'UserData', ida_user_data,...
                            'JacobianFn',@jacobianFn);
end

% Initialize IDA solver
% NOTE: Numeric values of t,x,xp will be input to diff_res(t,x,xp,ida_user_data) during IDA simulation. 
% Always use ida_user_data to define a function instead of defining a function inside diff_res.
IDAInit(@diff_res,t0,x0,xp0,options_IDA);

%% Find consistent ICs (only for IDACalcIC method)

if strcmp(param.init_Method,'IDACalcIC')
    try
        % Find consistent initial conditions
        [~, x0_initalized, ~] = IDACalcIC(t0+1e-3,'FindAlgebraic');
        x0 = x0_initalized';
    catch
        fprintf('Unsuccessful initialization: initialized algebraic states are NaN.\n')
        output.t = NaN;
        return
    end
elseif strcmp(param.init_Method,'SS')
    x0=x0';
end

%% Store ICs in output_sim structure
% -> For single-step method, the ICs will be the initial guess (OCV
% conditions) since consistent initial conditions are found during the 
% param.initime time frame once the simulation starts
% -> For IDACalcIC method, the consistent ICs have been determined in the
% section above. 

output_sim.t = t0;
output_sim.x = x0;
% Electrolyte concentration  
output_sim.ce = x0(param.ce_ind)*param.c0;    
% Positive concentration (based on SOC_init)
output_sim.csp = x0(param.csp_ind)*param.ctp;        
% Negative concentration (based on SOC_init)       
output_sim.csn = x0(param.csn_ind)*param.ctn;        
% Positive solid phase potential 
output_sim.phisp = x0(param.phisp_ind);         
% Negative solid phase potential
output_sim.phisn = x0(param.phisn_ind);         
% Liquid potential
output_sim.phie = x0(param.phie_ind);       
% Pore-wall flux0
output_sim.jp = x0(param.jp_ind);       
output_sim.jn = x0(param.jn_ind);      
% Current density
output_sim.curr_dens = x0(param.I_density_ind);
% Average and surface solid concentration
if strcmp(param.cs_discret,'FVM')
    [csp_surf, csp_avg, csp_avg_bulk, soc_p_bulk,...
        csn_surf, csn_avg, csn_avg_bulk, soc_n_bulk] =...
        calculateSolidConcentrationValues_FVM(output_sim.csp, output_sim.csn, output_sim.jp, output_sim.jn, param);
elseif strcmp(param.cs_discret,'FDM')
    [csp_surf, csp_avg, csp_avg_bulk, soc_p_bulk,...
        csn_surf, csn_avg, csn_avg_bulk, soc_n_bulk] =...
        calculateSolidConcentrationValues_FDM(output_sim.csp, output_sim.csn, param);
end
output_sim.csp_surf = csp_surf;
output_sim.csp_avg = csp_avg;
output_sim.csp_avg_bulk = csp_avg_bulk;
output_sim.soc_p_bulk = soc_p_bulk;
output_sim.csn_surf = csn_surf;
output_sim.csn_avg = csn_avg;
output_sim.csn_avg_bulk = csn_avg_bulk;
output_sim.soc_n_bulk = soc_n_bulk;
% Open-circuit voltage
output_sim.Up = param.Up(output_sim.csp_surf./param.ctp);
output_sim.Un = param.Un(output_sim.csn_surf./param.ctn);
% Overpotential (eta)
output_sim.etap = output_sim.phisp - output_sim.phie(1:Np) - output_sim.Up;
output_sim.etan = output_sim.phisn - output_sim.phie(Np+Ns+1:Np+Ns+Nn) - output_sim.Un;
% Voltage
switch param.VinterpMode
    case 1 % center of CV value 
        output_sim.voltage = output_sim.phisp(1) - output_sim.phisn(end) - param.Rc*(-output_sim.curr_dens)*param.Acell;
    case 2 % linear interpolation
        phis_pcc = 1.5*output_sim.phisp(1)-0.5*output_sim.phisp(2);
        phis_ncc = 1.5*output_sim.phisn(end)-0.5*output_sim.phisn(end-1);
        output_sim.voltage = phis_pcc - phis_ncc - param.Rc*(-output_sim.curr_dens)*param.Acell;
    case 3 % 3rd order hermite interpolation
        phis_pcc = electrodeBC_hermite_interp(output_sim.phisp,output_sim.curr_dens,param,'p');
        phis_ncc = electrodeBC_hermite_interp(output_sim.phisn,output_sim.curr_dens,param,'n');
        output_sim.voltage = phis_pcc - phis_ncc - param.Rc*(-output_sim.curr_dens)*param.Acell;
end

DAE_build_time = toc(tStart_DAE_build);
DAE_solve_time = 0;

%% Start IDA Solver

% Start at time t0
time = t0;

% If 'SS' include param.initime in tf
if strcmp(param.init_Method,'SS')
    tf = tf + param.initime;
end

% IDA solver while loop
while(time < tf)   
    % Check for simulation exit conditions
    [exit_indicator, exit_value, exit_message] =...
        simulationExitConditions(output_sim.voltage, output_sim.csp_surf, output_sim.csn_surf, param, output_sim.curr_dens);

    % Stop simulation since an exit condition has been reached
    if ~isempty(exit_indicator)
        if param.sim_print == 1
            fprintf('Simulation ending due to:\n');
            for i=1:length(exit_indicator)
                fprintf(exit_message{i});
            end
        end
        break
    end

    tStart_DAE_solve = tic;
    
    % IDA solver
    switch param.IDA_type 
        case 1
            try
                [~, time, x]   = IDASolve(tf,'OneStep');
            catch
                exit_indicator = {'IDASolve_error'};
                if param.sim_print == 1
                    fprintf('Simulation could not complete due to error in IDA solver\n');
                end
                break
            end
        case 2
            try
                time = time+param.deltaT;
                [~, time, x]   = IDASolve(time,'Normal');
            catch
                exit_indicator = {'IDASolve_error'};
                if param.sim_print == 1
                    fprintf('Simulation could not complete due to error in IDA solver\n');
                end
                break
            end
    end
    
    DAE_solve_time = DAE_solve_time+toc(tStart_DAE_solve);

    % Extract variables from x and store in output_sim
    output_sim = concatenateOutput(time, x', param, output_sim);

    if param.sim_print == 1
        fprintf(['t = ',num2str(time),' s\n']);
    end
end

% If no exit conditions triggered, make 'NaN' for string comparison
if isempty(exit_indicator)
    exit_indicator = {'NaN'};
end

% Store DAE build time and DAE solver time
output.DAE_build_time = DAE_build_time;
output.DAE_solve_time = DAE_solve_time;
output.exit_indicator = exit_indicator;
output.exit_value = exit_value;
output.exit_message= exit_message;

%% Keep original simulation (this may make output structure to have larger file size-suppress to save memory)
 % For the single-step method case, output.original will be the solution containing the initialization part
if param.save_memory == 1
    output.original = output_sim;
end

%%  Remove the initialization solution (only for single-step method)

if strcmp(param.init_Method,'SS')
    % If the initialization finished (no exit conditions triggered during initialization)
    if time >= param.initime
        % Find index of output without initialization part 
        pos_ind = find(output_sim.t-param.initime>=0);
        % Include index right before time turns positive so that we can interpolate to find t=0 later
        pos_ind = [pos_ind(1)-1; pos_ind];
        % Store the solution without initialization part in output_sim
        output_sim.t = output_sim.t(pos_ind)-param.initime;
        output_sim.x = output_sim.x(pos_ind,:);
        output_sim.ce = output_sim.ce(pos_ind,:);
        output_sim.csp = output_sim.csp(pos_ind,:);
        output_sim.csn = output_sim.csn(pos_ind,:);
        output_sim.phisp = output_sim.phisp(pos_ind,:);
        output_sim.phisn = output_sim.phisn(pos_ind,:);
        output_sim.phie = output_sim.phie(pos_ind,:);
        output_sim.jp = output_sim.jp(pos_ind,:);
        output_sim.jn = output_sim.jn(pos_ind,:);
        output_sim.curr_dens = output_sim.curr_dens(pos_ind);
        output_sim.csp_surf = output_sim.csp_surf(pos_ind,:);
        output_sim.csp_avg = output_sim.csp_avg(pos_ind,:);
        output_sim.csp_avg_bulk = output_sim.csp_avg_bulk(pos_ind,:);
        output_sim.soc_p_bulk = output_sim.soc_p_bulk(pos_ind);
        output_sim.csn_surf = output_sim.csn_surf(pos_ind,:);
        output_sim.csn_avg = output_sim.csn_avg(pos_ind,:);
        output_sim.csn_avg_bulk = output_sim.csn_avg_bulk(pos_ind,:);
        output_sim.soc_n_bulk = output_sim.soc_n_bulk(pos_ind);
        output_sim.voltage = output_sim.voltage(pos_ind);
        output_sim.Up = output_sim.Up(pos_ind,:);
        output_sim.Un = output_sim.Un(pos_ind,:);
        output_sim.etap = output_sim.etap(pos_ind,:);
        output_sim.etan = output_sim.etan(pos_ind,:);
    % If the initialization did not finish (exit conditions triggered during initialization)
    % -> then make the solution that triggered the exit condition be the solution at t=0
    else
        output_sim.t = 0;
        output_sim.x = output_sim.x(end,:);
        output_sim.ce = output_sim.ce(end,:);
        output_sim.csp = output_sim.csp(end,:);
        output_sim.csn = output_sim.csn(end,:);
        output_sim.phisp = output_sim.phisp(end,:);
        output_sim.phisn = output_sim.phisn(end,:);
        output_sim.phie = output_sim.phie(end,:);
        output_sim.jp = output_sim.jp(end,:);
        output_sim.jn = output_sim.jn(end,:);
        output_sim.curr_dens = output_sim.curr_dens(end);
        output_sim.csp_surf = output_sim.csp_surf(end,:);
        output_sim.csp_avg = output_sim.csp_avg(end,:);
        output_sim.csp_avg_bulk = output_sim.csp_avg_bulk(end,:);
        output_sim.soc_p_bulk = output_sim.soc_p_bulk(end);
        output_sim.csn_surf = output_sim.csn_surf(end,:);
        output_sim.csn_avg = output_sim.csn_avg(end,:);
        output_sim.csn_avg_bulk = output_sim.csn_avg_bulk(end,:);
        output_sim.soc_n_bulk = output_sim.soc_n_bulk(end);
        output_sim.voltage = output_sim.voltage(end);
        output_sim.Up = output_sim.Up(end,:);
        output_sim.Un = output_sim.Un(end,:);
        output_sim.etap = output_sim.etap(end,:);
        output_sim.etan = output_sim.etan(end,:);
    end
end

%% Solution interpolation/clean up

% If 'SS' remove param.initime in tf 
if strcmp(param.init_Method,'SS')
    tf = tf - param.initime;
end

if length(output_sim.t)>1
    % Find actual final time if exit triggered by tf or voltage cutoffs
    % This is because the IDA solver will output the solution one time step
    % after the exit condition is triggered. We need to manually interpolate
    % the last solution to find the last solution corresponding to the trigger condition.
    
    % Lower cutoff voltage triggered
    if strcmp(exit_indicator{1},'V_cutoffLower') 
        time_cutoff = interp1(output_sim.voltage,output_sim.t,param.lowerCutoffVoltage);
    % Upper cutoff voltage triggered
    elseif strcmp(exit_indicator{1},'V_cutoffUpper') 
        time_cutoff = interp1(output_sim.voltage,output_sim.t,param.upperCutoffVoltage);
    % cs_surf > cs_max (square root in BV equation will be negative)
    elseif strcmp(exit_indicator{1},'cspSurf_aboveCspMax') || strcmp(exit_indicator{1},'csnSurf_aboveCsnMax')   
        time_cutoff = output_sim.t(end);
    % Simulated ended prematurely due to issues in IDA solver
    elseif strcmp(exit_indicator,'IDASolve_error') 
        time_cutoff = output_sim.t(end);
    % Exit triggered by simulation time being greater than tf
    elseif time >= tf % same if condition as strcmp(exit_indicator{1},'NaN') 
        time_cutoff = tf;
    end

    % Store initial states data for future cycles
    output.x_initials = interp1(output_sim.t,output_sim.x,time_cutoff);

    % deltaT for interpolation
    deltaT = param.deltaT_sim; 
    % Interpolation time
    time_interp = (t0:deltaT:time_cutoff)';
    % Interpolate solution so that we have finer deltaT points
    output.t = time_interp;
    ce_interp = interp1(output_sim.t, output_sim.ce, time_interp);
    csp_interp = interp1(output_sim.t, output_sim.csp, time_interp);
    csn_interp = interp1(output_sim.t, output_sim.csn, time_interp);
    phisp_interp = interp1(output_sim.t, output_sim.phisp, time_interp);
    phisn_interp = interp1(output_sim.t, output_sim.phisn, time_interp);
    phie_interp = interp1(output_sim.t, output_sim.phie, time_interp);
    jp_interp = interp1(output_sim.t, output_sim.jp, time_interp);
    jn_interp = interp1(output_sim.t, output_sim.jn, time_interp);
    curr_dens_interp = interp1(output_sim.t, output_sim.curr_dens, time_interp);
    csp_surf_interp = interp1(output_sim.t, output_sim.csp_surf, time_interp);
    csp_avg_interp = interp1(output_sim.t, output_sim.csp_avg, time_interp);
    csp_avg_bulk_interp = interp1(output_sim.t, output_sim.csp_avg_bulk, time_interp);
    soc_p_bulk_interp = interp1(output_sim.t, output_sim.soc_p_bulk, time_interp);
    csn_surf_interp = interp1(output_sim.t, output_sim.csn_surf, time_interp);
    csn_avg_interp = interp1(output_sim.t, output_sim.csn_avg, time_interp);
    csn_avg_bulk_interp = interp1(output_sim.t, output_sim.csn_avg_bulk, time_interp);
    soc_n_bulk_interp = interp1(output_sim.t, output_sim.soc_n_bulk, time_interp);
    voltage_interp = interp1(output_sim.t, output_sim.voltage, time_interp);
    Up_interp = interp1(output_sim.t, output_sim.Up, time_interp);
    Un_interp = interp1(output_sim.t, output_sim.Un, time_interp);
    etap_interp = interp1(output_sim.t, output_sim.etap, time_interp);
    etan_interp = interp1(output_sim.t, output_sim.etan, time_interp);
% This is for a case when the output solution has only one entry
% -> this means that an exit condition was triggered at the start of the simulation
else 
    output.t = output_sim.t;
    ce_interp = output_sim.ce;
    csp_interp = output_sim.csp;
    csn_interp = output_sim.csn;
    phisp_interp = output_sim.phisp;
    phisn_interp = output_sim.phisn;
    phie_interp = output_sim.phie;
    jp_interp = output_sim.jp;
    jn_interp = output_sim.jn;
    curr_dens_interp = output_sim.curr_dens;
    csp_surf_interp = output_sim.csp_surf;
    csp_avg_interp = output_sim.csp_avg;
    csp_avg_bulk_interp = output_sim.csp_avg_bulk;
    soc_p_bulk_interp = output_sim.soc_p_bulk;
    csn_surf_interp = output_sim.csn_surf;
    csn_avg_interp = output_sim.csn_avg;
    csn_avg_bulk_interp = output_sim.csn_avg_bulk;
    soc_n_bulk_interp = output_sim.soc_n_bulk;
    voltage_interp = output_sim.voltage;
    Up_interp = output_sim.Up;
    Un_interp = output_sim.Un;
    etap_interp = output_sim.etap;
    etan_interp = output_sim.etan;

    % Store initial states data for future cycles
    output.x_initials = output_sim.x;
end

% Define output structure
output.ce = ce_interp;
output.csp = csp_interp;
output.csn = csn_interp;
output.phisp = phisp_interp;
output.phisn = phisn_interp;
output.phie = phie_interp;
output.jp = jp_interp;
output.jn = jn_interp;
output.curr_dens = curr_dens_interp;
output.csp_surf = csp_surf_interp;
output.csp_avg = csp_avg_interp;
output.csp_avg_bulk= csp_avg_bulk_interp;
output.soc_p_bulk = soc_p_bulk_interp;
output.csn_surf = csn_surf_interp;
output.csn_avg = csn_avg_interp;
output.csn_avg_bulk = csn_avg_bulk_interp;
output.soc_n_bulk = soc_n_bulk_interp;
output.voltage = voltage_interp;
output.Up = Up_interp;
output.Un = Un_interp;
output.etap = etap_interp;
output.etan = etan_interp;
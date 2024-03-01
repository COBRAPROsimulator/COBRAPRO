---
title: 'COBRAPRO: A MATLAB toolbox for Physics-based Battery Modelling and Co-simulation Parameter Optimization'
tags:
  - lithium-ion batteries
  - physics-based modeling
  - DFN model
  - P2D model
  - parameter identification
  - particle swarm optimization
authors:
  - name: Sara Ha
    orcid: 0009-0005-9878-3537
    affiliation: 1 
  - name: Simona Onori
    orcid: 0000-0002-6556-2608
    affiliation: 2
affiliations:
 - name: Mechanical Engineering, Stanford University, 440 Escondido Mall, Stanford 94305, CA, USA
   index: 1
 - name: Energy Science and Engineering, Stanford University, 367 Panama Mall, Stanford 94305, CA, USA
   index: 2
date: 01 March 2024
bibliography: paper.bib
---

# Summary
COBRAPRO (**Co**-simulation **B**atte**r**y Modeling for **A**ccelerated **P**a**r**ameter **O**ptimization) is a physics-based battery modeling software with the capability to perform closed-loop parameter optimization using experimental data. COBRAPRO is based on the Doyle-Fuller-Newman (DFN) model [1], which is most widely-accepted high-fidelity model that considers the lithium-ion transport and charge conservation in the liquid electrolyte and solid electrodes, and kinetics at the solid and liquid interface during lithium intercalation and deintercalation. Such physics-based models have found applications in battery design [2], [3] and advanced battery management systems to ensure reliable and safe operation of electric vehicles [4]. The DFN model encompasses several physical parameters, such as geometric, stoichiometric, concentration, transport, and kinetic parameters, which are often unknown and need to be determined to accurately predict battery response under various usage scenarios. Direct measurement through cell tear-down experiments is a viable but labor-intensive process [5], [6], [7]. Furthermore, parameters obtained through experimental characterization may not be suitable for the DFN model [7], as the model is a simplified representation of a real battery, assuming perfectly spherical particles, neglecting electrode heterogeneity, and only considering internal dynamics in one dimension. With COBRAPRO, users can noninvasively identify parameters for any given battery using readily available current-voltage data from a battery cycler. COBRAPRO optimizes the DFN parameters by minimizing the error between the simulated and experimentally observed data through an embedded parameter optimization routine.

# Statement of need

Currently, battery DFN modeling tools lack integrated identification routines and primarily rely on parameter values obtained from literature [8] [9] [11] [10]. COMSOL Multiphysics [cite], a commercially available finite element modeling software utilized for DFN simulations, features a Battery Design Module with a parameter library for various chemistries. However, inherent material inconsistencies and manufacturing cell-to-cell variability may render parameters provided by COMSOL unsuitable for modeling a battery even with the same chemistry and electrolyte material. Additionally, COMSOL lacks a built-in parameter identification feature. To address this, one can use COMSOL’s LiveLink for MATLAB to establish communication between the two software platforms such that results from COMSOL can be transferred to MATLAB to leverage the versatile suite of optimizers in MATLAB [12]. Nonetheless, the expensive licensing fee and proprietary nature of COMSOL create barriers to public access, limiting collaboration and reproducibility.

In contrast, several open-source DFN model simulation tools have been released such as PyBaMM [8], LIONSIMBA [9], PETLION [13], DEARLIBS [14], fastDFN, and MPET [10]. Among these packages, DEARLIBS is the only software with the capability to perform closed-loop parameter identification. Other packages rely on literature-derived parameter values to generate simulation results. Taking inspiration from DEARLIBS, COBRAPRO aims to address three primary challenges in the DFN model.

## Challenge 1. Computational complexity 
- **Issue:** The DFN model is also called the pseudo-two-dimensional (P2D) due to the coupling of the cell thickness (x-direction) and radial particle (r-direction) dimensions, contributing to the computational complexity of the system.
- **Solution:** COBRAPRO uses a fast solver that significantly improves the model computation speed compared to DEARLIBS. For 10 discretized points in each domain of the cell (cathode, separator, anode, and positive and negative active material particles) at 1C discharge, COBRAPRO solves the DFN model in 0.708 seconds, while DEARLIBS takes 2.54 minutes, which is a two orders of magnitude improvement (~257 times). LIONSIMBA took 1.13 seconds and PyBaMM took 0.237 seconds, which are comparable to the COBRAPRO’s computation time. For larger discretization points, we observed up to three orders of magnitude improvement in computation speed from COBRAPRO to DEARLIBS.

## Challenge 2. Consistent initial conditions
- **Issue:** The partial differential equations (PDEs) governing the DFN model are discretized in the x and r directions to form a system of ordinary differential equations (ODEs) and algebraic equations (AEs), also called differential-algebraic equations (DAEs). To solve the DAE system, the correct AEs are required, which are typically not known a priori for the DFN model. Inconsistent initial conditions result in either a failure to start the simulation or the model diverging towards an incorrect solution [15].
- **Solution:** The single step approach [16], a robust initialization method, is implemented in COBRAPRO that automatically determines the initial conditions and seamlessly simulates the DFN model.

## Challenge 3. Unknown model parameters
- **Issue:** As highlighted earlier, battery parameters are frequently unknown, and even if obtained through experimental characterization, parameter calibration is essential to accurately model the battery.
- **Solution:** A co-simulation parameter optimization framework is developed that determines the parameters by minimizing the cost function, defined in terms of the error between the experimental and simulated voltage and state-of-charge curves. The particle swarm optimization (PSO), a gradient-free population-based algorithm, is employed due to its suitability for nonlinear models like the DFN model. COBRAPRO employs MATLAB’s Parallel Computing Toolbox, accelerating PSO through multicore processing.

# Core Capabilities
-	DFN model implementation using finite volume method (FVM) discretization and SUNDIALS IDA solver
-	Radial solid particle discretization options:
  - FVM (3rd order Hermite interpolation is utilized to accurately estimate the particle surface concentration to account for the sharp concentration gradients near the particle surface)
  -	Finite difference method (FDM) 
- DAE initialization options:
  -	Single-step approach [16]
  -	SUNDIALS IDACalcIC
-	Parameter identification routine using experimental current-voltage data (Examples/Parameter_Identification_Routines)
-	Simulating battery cycling 
  - Constant current cycling 
  - Hybrid pulse power characterization (HPPC) test simulation 
  - Dynamic current simulation 
-	Local sensitivity analysis 

Visit COBRAPRO’s Github page for example codes that demonstrate each of the features listed above.

# References

---
title: 'COBRAPRO: A MATLAB toolbox for Physics-based Battery Modeling and Co-simulation Parameter Optimization'
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
COBRAPRO (**Co**-simulation **B**atte**r**y Modeling for **A**ccelerated **P**a**r**ameter **O**ptimization) is a physics-based battery modeling software with the capability to perform closed-loop parameter optimization using experimental data. COBRAPRO is based on the Doyle-Fuller-Newman (DFN) model [@doyle_modeling_1993], which is most widely-accepted high-fidelity model that considers the lithium-ion transport and charge conservation in the liquid electrolyte and solid electrodes, and kinetics at the solid and liquid interface during lithium intercalation and deintercalation. Such physics-based models have found applications in battery design [@dai_graded_2016], [@couto_lithiumion_2023] and advanced battery management systems to ensure reliable and safe operation of electric vehicles [@kolluri_realtime_2020]. The DFN model encompasses several physical parameters, such as geometric, stoichiometric, concentration, transport, and kinetic parameters, which are often unknown and need to be determined to accurately predict battery response under various usage scenarios. Direct measurement through cell tear-down experiments is a viable but labor-intensive process [@ecker_parameterization_2015a], [@schmalstieg_full_2018a], [@chen_development_2020]. Furthermore, parameters obtained through experimental characterization may not be suitable for the DFN model [@chen_development_2020]. The DFN model simplifies the representation of a real battery by assuming perfectly spherical particles, neglecting electrode heterogeneity, and considering internal dynamics in only one dimension. With COBRAPRO, users can noninvasively identify parameters for any given battery using readily available current-voltage data from a battery cycler. COBRAPRO optimizes the DFN parameters by minimizing the error between the simulated and experimentally observed data through an embedded parameter optimization routine.

# Statement of need

Currently, battery DFN modeling tools lack integrated identification routines and primarily rely on parameter values obtained from literature [@sulzer_python_2021] [@torchio_lionsimba_2016] [@smith_multiphase_2017] [@berliner_methods_2021]. COMSOL Multiphysics&copy; [@comsol], a commercially available finite element modeling software utilized for DFN simulations, features a Battery Design Module [@comsol_battery] with a parameter library for various chemistries. However, inherent material inconsistencies and manufacturing cell-to-cell variability [@harris_failure_2017] may render parameters provided by COMSOL unsuitable for modeling a battery even with the same chemistry and electrolyte material. Additionally, COMSOL lacks a built-in parameter identification feature. To address this, one can use COMSOL’s *LiveLink&trade; for MATLAB&copy;* to establish communication between the two software platforms such that results from COMSOL can be transferred to MATLAB to leverage the versatile suite of optimizers in MATLAB [@pozzato_general_2023]. Nonetheless, the expensive licensing fee and proprietary nature of COMSOL create barriers to public access, limiting collaboration and reproducibility.

In contrast, several open-source DFN model simulation tools have been released such as PyBaMM [@sulzer_python_2021], LIONSIMBA [@torchio_lionsimba_2016], PETLION [@berliner_methods_2021], DEARLIBS [@lee_robust_2021], fastDFN[@fastDFN], and MPET [@smith_multiphase_2017]. Among these packages, DEARLIBS is the only software with the capability to perform closed-loop parameter identification using experimental data. Other packages rely on literature-derived parameter values to generate simulation results. Taking inspiration from DEARLIBS, COBRAPRO aims to address three primary challenges in the DFN model.

## Challenge 1. Computational complexity 
- **Issue:** The DFN model is also called the pseudo-two-dimensional (P2D) due to the coupling of the cell thickness (x-direction) and radial particle (r-direction) dimensions, contributing to the computational complexity of the system.
- **Solution:** COBRAPRO uses a fast solver that significantly improves the model computation speed compared to DEARLIBS. For 10 discretized points in each domain of the cell (positive and negative electrodes, separator, and positive and negative active material particles) at 1C discharge, COBRAPRO solves the DFN model in 0.708 seconds, while DEARLIBS takes 2.54 minutes, which is a two orders of magnitude improvement (~257 times). Under the same simulation conditions, LIONSIMBA and PyBaMM computed the model in 1.13 seconds and 0.237 seconds, respectively, which are comparable to the COBRAPRO’s computation time. For larger discretization points, we observed up to three orders of magnitude improvement in computation speed from COBRAPRO to DEARLIBS.

## Challenge 2. Consistent initial conditions
- **Issue:** The partial differential equations (PDEs) governing the DFN model are discretized in the x- and r-directions to form a system of ordinary differential equations (ODEs) and algebraic equations (AEs), also called differential-algebraic equations (DAEs). To solve the DAE system, the correct AEs are required, which are typically not known *a priori* for the DFN model. Inconsistent initial conditions result in either a failure to start the simulation or the model diverging towards an incorrect solution [@methekar_perturbation_2011].
- **Solution:** The single step approach [@lawder_extending_2015], a robust initialization method, is implemented in COBRAPRO that automatically determines the initial conditions and seamlessly simulates the DFN model.

## Challenge 3. Unknown model parameters
- **Issue:** As highlighted earlier, battery parameters are frequently unknown, and even if obtained through experimental characterization, parameter calibration is essential to accurately model the battery.
- **Solution:** A co-simulation parameter optimization framework is developed that determines the parameters by minimizing the cost function, defined in terms of the error between the experimental and simulated voltage and state-of-charge curves. The particle swarm optimization (PSO), a gradient-free population-based algorithm, is employed due to its suitability for nonlinear models like the DFN model. COBRAPRO employs MATLAB’s Parallel Computing Toolbox, accelerating PSO through multicore processing.

# Core Capabilities
- **Parameter identification routine:** Utilizes PSO to optimize parameters using experimental current-voltage data
- **DFN model implementation:** Finite volume method (FVM) discretization of the PDEs to form a DAE and SUNDIALS IDA solver to solve the DAE system
- **Solid particle radial discretization options:**
  - FVM (3rd order Hermite interpolation is utilized to accurately estimate the particle surface concentration to account for the sharp concentration gradients near the particle surface [@xu_comparative_2023])
  - Finite difference method (FDM)
- **DAE initialization options:**
  - Single-step approach [@lawder_extending_2015]
  - SUNDIALS IDACalcIC
- **Simulating battery cycling:**
  - Constant current profiles
  - Hybrid pulse power characterization (HPPC) profiles
  - Dynamic current profiles 
- **Local sensitivity analysis:** Perturbs parameters around specific reference values to determine sensitive parameters for a given current profile

# Examples 
To get started, view example codes included in COBRAPRO's ```Examples``` folder.
- ``Examples/Cycling``: examples showing how to perform battery cycling simulations using experimentally identified parameters


# Acknowledgements
The authors thank the Bits and Watts Initiative within the Precourt Institute for Energy at Stanford University for its partial financial support. We thank Dr. Le Xu for all the insightful discussions that greatly contributed to the enhancement of COBRAPRO. We extend our thanks to Alexis Geslin, Joseph Lucero, and Maitri Uppaluri for testing COBRAPRO and providing valuable feedback.

# References

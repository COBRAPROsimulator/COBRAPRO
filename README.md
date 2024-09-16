# COBRAPRO: A MATLAB toolbox for Physics-based Battery Modeling and Co-simulation Parameter Optimization

<a name="readme-top"></a>

<!--[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]-->

COBRAPRO (Co-simulation Battery Modeling for Accelerated Parameter Optimization) is a MATLAB software for physics-based modeling of lithium-ion batteries (LIB) with an embedded parameter identification routine. We aim to provide the battery modeling community with a versatile toolbox for calibrating battery models, a crucial process to achieve accurate simulation results for predicting real-world battery responses under various operating conditions. 

Please refer to the following publication for more details on the numerical methods, determination of consistent initial condition, and parameter optimization pipeline proposed in CORBAPRO:

S. Ha and S. Onori, “COBRAPRO: An Open-Source Software for the Doyle-Fuller-Newman Model with Co-Simulation Parameter Optimization Framework,” J. Electrochem. Soc., Aug. 2024, doi: 10.1149/1945-7111/ad7292.

## Table of contents ##

  * [What is COBRAPRO?](#toc1)
  * [Why COBRAPRO?](#toc2)
  * [System requirements](#toc3)
  * [Installation](#toc4)
  * [Testing](#toc5)
  * [Examples](#toc6)
  * [Contributing](#toc7)
  * [Known issues](#toc8)

## What is COBRAPRO? <a name="toc1"></a> ##

COBRAPRO implements the Dolye-Fuller-Newman (DFN) model, also known as the pseudo-two-dimensional (P2D) model, which is a high-fidelity LIB model considering the lithium-ion mass and charge conservation in the liquid electrolyte and solid electrodes, and Butler-Volmer kinetics. Parameter calibration, or identification, is a primary challenge in implementing the DFN model since the parameters such as geometric, transport, kinetic, concentration, and stoichiometric are often not known _a prioi_. 

In response to this challenge, COBRAPRO allows users to identify parameters of any battery cells based on their experimental current-voltage profiles. COBRAPRO solves an optimization problem that minimizes the error between the experimental and simulated voltage and state-of-charge curves to identify the parameters of interest. Although the software employs particle swarm optimization (PSO) by default, users have the flexibility modify the code to implement other MATLAB optimization algorithms such as `ga`, `fmincon`, `patternsearch`, and more. 

## Why COBRAPRO? <a name="toc2"></a> ##

Compared to currently available DFN open-source packages such as PyBaMM, DEARLIBS, LIONSIMBA, PETION, fastDFN, and MPET, DEARLIBS and COBRAPRO are the only codes with an integrated identification routine. Given the need for numerous model simulations during parameter optimization, achieving efficient computation time is critical. COBRAPRO addresses this need with a fast solver and PSO parallel computing, resulting in model simulations up to three orders of magnitude faster than DEARLIBS and accelerated PSO through multicore processing.

## Software dependencies <a name="toc3"></a> ##
* MATLAB 2018b and later
* MATLAB Global Optimization Toolbox
*	MATLAB Parallel Computing Toolbox
*	SUNDIALS 2.6.2
*	CasADi (MATLAB version)
*	Xcode (for macOS users only)
*	MinGW (for Window users only)

Installation section below shows how to install the required software.

## Installation <a name="toc4"></a> ##

1. Download COBRAPRO by downloading the zip file or cloning this repository by typing in Terminal:
   ```
   git clone https://github.com/COBRAPROsimulator/COBRAPRO.git
   ```
2. Download [MATLAB](https://www.mathworks.com/downloads/) if not installed already. Make sure to select both the Global Optimization Toolbox and Parallel Computing Toolbox during the installation process.
   
3. Download [SUNDIALS 2.6.2](https://computing.llnl.gov/sites/default/files/inline-files/sundials-2.6.2.tar.gz) and unzip the folder. Relocate the sundials-2.6.2 folder inside the COBRAPRO folder.

4. Download the latest version of [CasADi (MATLAB version)](https://web.casadi.org/get/) corresponding to your operating system. Unzip and move your CasADi folder inside the COBRAPRO folder. Your COBRAPRO folder should now contain the sundials-2.6.2 and CasADi folders.

5. Before we can install SUNDIALS, the following software are required to compile the mex files that will interface with the SUNDIALS IDA solver:
   - __Mac users__: Download [Xcode](https://developer.apple.com/xcode/) application (can be downloaded from Apple’s App Store). Once Xcode[^1] is installed, proceed to accept the license agreement. This can be done by opening the Xcode application, which will launch a license agreement window and click the “Agree” icon, or type
     ```
     sudo xcodebuild -license accept
     ```
     in Terminal. If the license is not accepted, MATLAB may give an error such as “Xcode is installed, but its license has not been accepted”.
   - __Window users__: Download [MinGW](https://www.mathworks.com/matlabcentral/fileexchange/52848-matlab-support-for-mingw-w64-c-c-fortran-compiler)
     
6. Now you are ready to run `install_COBRAPRO.m`, which is located in the main COBRAPRO folder. `install_COBRAPRO.m` will install SUNDIALS by calling the `sundials-2.6.2/sundialsTB/install_STB.m` file and automatically add the required folders to your MATLAB path. Run `install_COBRAPRO.m` and respond to the prompts displayed in the Command Window in the following manner:

```
MEX files will be compiled and built using the above options
   Proceed? (y/n)
```
&rarr; Type `y` and hit enter

```
Compile CVODES interface? (y/n)
```
&rarr; Type `n` and hit enter
```
Compile IDAS  interface? (y/n)
```
&rarr; Type `y` and hit enter
```
Compile KINSOL  interface? (y/n)
```
&rarr; Type `n` and hit enter
```
MEX files were successfully created.
    Install toolbox? (y/n) 
```
&rarr; Type `y` and hit enter
```
Specify the location where you wish to install the toolbox.
The toolbox will be installed in a subdirectory "sundialsTB".
Enter return to cancel the installation.
Installation directory:
```
&rarr; Just hit enter
```
Type the name of your CasADi folder (case-sensitive):
```
&rarr; Type the name of the CasADi folder exactly as it appears and hit enter

7. Successful installation will output to the Command Window:
```
COBRAPRO installed successfully! To get started, try running scripts in the Examples folder.
```
## Testing <a name="toc5"></a> ##

Automated test codes are provided in `test` folder:
- `test_1_casadiCheck.m` checks that CasADi is installed and working properly. Successful run will output to Command Window:
```
test_1 successful: CasADi is working properly!
```
- `test_2_comsolValidation.m` validates COBRAPRO against results generated from COMSOL Multiphysics[^2] as a benchmark. This ensures that COBRAPRO is installed properly and that the SUNDIALS IDA solver is working as expected. Successful validation will output to Command Window:
```
test_2 successful: COBRAPRO is working as expected! Results validated against COMSOL.
```
- `test_3_psoCheck.m` ensures that MATLAB's Global Optimization Toolbox and Parallel Computing Toolbox are installed, and makes sure that the PSO in parallel (required for parameter identification) is working correctly. Successful run will output to Command Window:
```
test_3 successful: PSO in parallel is working as expected!
```

## Examples <a name="toc6"></a> ##
In the ```Examples``` folder, you will find example codes that will help you get started.
* ```Examples/Cycling```: examples showing how to perform battery cycling simulations using experimentally identified parameters 
  * ```Examples/Cycling/cycle_CC.m```: simulating constant current (CC) cycling experiments and result visualization (voltage, state-of-charge, internal variable curves)
  * ```Examples/Cycling/cycle_HPPC.m```: simulating hybrid pulse power characterization (HPPC) profile and result visualization (voltage, state-of-charge, internal variable curves)
  * ```Examples/Cycling/cycle_UDDS.m```: simulating driving cycle profile and result visualization (voltage, state-of-charge, internal variable curves)
* ```Examples/Parameter_Identification_Routines```: examples showing how to perform parameter identification using PSO
&rarr; NOTE: In general, these scripts take a while to run. Using a processor with multiple cores, e.g., 12 or 24 cores, will significantly speed up the PSO. Also, PSO particle size and PSO exit conditions affect the PSO convergence accuracy and time.)
  * ```Examples/Parameter_Identification_Routines/DFN_pso_0_05C.m```: parameter identification using C/20 discharge data
  * ```Examples/Parameter_Identification_Routines/DFN_pso_HPPC.m```: parameter identification using HPPC data (given same number of PSO particles and PSO exit conditions, takes longer to run than ```DFN_pso_0_05C.m``` since HPPC takes much longer to run than C/20 discharge)
* ```Examples/Parameter_Identification_Results```: examples showing parameter identification results
  * ```Examples/Parameter_Identification_Results/DFN_pso_0_05C_identification.m```: parameter identification results using C/20 discharge data 
  * ```Examples/Parameter_Identification_Results/DFN_pso_HPPC_identification.m```: parameter identification results using HPPC data
  * ```Examples/Parameter_Identification_Results/DFN_pso_UDDS_validation.m```: parameter identification validation using UDDS data
* ```Examples/Local_Sensitivity_Analysis```: examples showing how to perform local sensitivity analysis (LSA)
  * ```Examples/Local_Sensitivity_Analysis/DFN_LSA_Corr_CC.m```: LSA and correlation analysis on CC profile
  * ```Examples/Local_Sensitivity_Analysis/DFN_LSA_Corr_HPPC.m```: LSA and correlation analysis on HPPC profile

## Contributing <a name="toc7"></a> ##
We welcome contributions from the community to improve COBRAPRO!
* To report bugs, ask questions, and get help, please open a new issue through the Github issues page. Be as specific as possible (including screenshots, sample codes) for efficient communication. 
* To make changes to the code or add new functions, 1) fork the repo and create your branch from main, 2) make your changes to the code, and 3) open a Pull request. Once approved, your contribution will be merged into the master branch.
* For general discussions and project ideas, open a new Discussions through the Github issues page. You can also contact Sara Ha (<sungyeon.sara.ha@stanford.edu>).
  
## Known issues <a name="toc8"></a> ##
1. To run COBRAPRO, only the SUNDIALS IDAS interface is required. In Installation step 4, if you install the KINSOL interface, you may run into the following issue:
```
Error using mex
COBRAPRO/sundials-2.6.2/sundialsTB/kinsol/kim/src/kim.c:687:24: error: non-void function 'KIM_Stats' should return a value [-Wreturn-type]
if (kimData == NULL) return;
COBRAPRO/sundials-2.6.2/sundialsTB/kinsol/kim/src/kim.c:687:24: error: non-void function 'KIM_Free' should return a value [-Wreturn-type]
return;
2 errors generated.
```
To fix this issue, please go to `sundials-2.6.2/sundialsTB/kinsol/kim/src/kim.c` and modify line 687 to
```
if (kimData == NULL) return NULL;
```
and modify line 815 to
```
return NULL;
```
2. If you happen to run into the following error when installing/running COBRAPRO:
```
Error using mex
'idm.mexmaca64' locked by mexLock API.
``` 
This seems to occur when MATLAB is trying to compile the mex files and mexLock is triggered. To resolve the issue, please restart MATLAB. 

[^1]: Note that Xcode requires ~3.4 GB of storage space.  
[^2]: COMSOL Multiphsyics is a commerically available finite element analysis software.

## How to cite this code
If you use this code in your research, please cite our JOSS paper and the accompanying JES paper:
```
@article{cobrapro_joss_2024,
   author = {Ha, Sara and Onori, Simona},
   doi = {arXiv:2404.10022},
   journal = {Journal of Open Source Software},
   year = {2024}
   title = {{COBRAPRO: A MATLAB toolbox for Physics-based Battery Modeling and Co-simulation Parameter Optimization}},
}
```
```
@article{cobrapro_jes_2024,
   author = {Ha, Sara and Onori, Simona},
   doi = {10.1149/1945-7111/ad7292},
   journal = {Journal of The Electrochemical Society},
   year = 2024,
   title = {{COBRAPRO: An Open-Source Software for the Doyle-Fuller-Newman Model with Co-Simulation Parameter Optimization Framework}},
} 
```

## Contributors

[![All Contributors](https://img.shields.io/github/all-contributors/COBRAPROsimulator/COBRAPRO?color=ee8449&style=flat-square)](#contributors)

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->

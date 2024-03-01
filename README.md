# COBRAPRO: A MATLAB toolbox for Physics-based Battery Modelling and Co-simulation Parameter Optimization

<a name="readme-top"></a>

<!--[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]-->

COBRAPRO (Co-simulation Battery Modeling for Accelerated Parameter Optimization) is a MATLAB software for physics-based modeling of lithium-ion batteries (LIB) with an embedded parameter identification routine. We aim to provide the battery modeling community with a versatile toolbox for calibrating battery models, a crucial process to achieve accurate simulation results for predicting real-world battery responses under various operating conditions.

## Table of contents ##

  * [What is COBRAPRO?](#toc1)
  * [Why COBRAPRO?](#toc2)
  * [System requirements](#toc3)
  * [Installation](#toc4)
  * [Testing](#toc5)
  * [Examples](#toc6)
  * [Contributing](#toc7)
  * [Known issues during installation](#toc8)

## What is COBRAPRO? <a name="toc1"></a> ##

COBRAPRO implements the Dolye-Fuller-Newman (DFN) model, also known as the pseudo-two-dimensional (P2D) model, which is a high-fidelity LIB model considering the lithium-ion mass and charge conservation in the liquid electrolyte and solid electrodes, and Butler-Volmer kinetics. Parameter calibration, or identification, is a primary challenge in implementing the DFN model since the parameters such as geometric, transport, kinetic, concentration, and stoichiometric are often not known _a prioi_. 

In response to this challenge, COBRAPRO allows users to identify parameters of any battery cells based on their experimental current-voltage profiles. COBRAPRO solves an optimization problem that minimizes the error between the experimental and simulated voltage and state-of-charge curves to identify the parameters of interest. Although the software employs particle swarm optimization (PSO) by default, users have the flexibility modify the code to implement other MATLAB optimization algorithms such as `ga`, `fmincon`, `patternsearch`, and more. 

## Why COBRAPRO? <a name="toc2"></a> ##

Compared to currently available DFN open-source packages such as PyBaMM, DEARLIBS, LIONSIMBA, PETION, fastDFN, and MPET, DEARLIBS and COBRAPRO are the only codes with an integrated identification routine. Given the need for numerous model simulations during parameter optimization, achieving efficient computation time is critical. COBRAPRO addresses this need with a fast solver and PSO parallel computing, resulting in model simulations up to three orders of magnitude faster than DEARLIBS and accelerated PSO through multicore processing.

## Software dependencies <a name="toc3"></a> ##
* MATLAB 2018b and later
*	MATLAB Global Optimization Toolbox
*	MATLAB Parallel Computing Toolbox
*	SUNDIALS 2.6.2
*	CasADi
*	Xcode (for macOS users only)
*	MinGW (for Window users only)

Installation section below shows how to install the required software.

## Installation <a name="toc4"></a> ##

1. Download COBRAPRO by downloading the zip file or cloning this repository by typing in Terminal:
   ```
   git clone https://github.com/COBRAPROsimulator/COBRAPRO.git
   ```
2. Download [SUNDIALS 2.6.2](https://computing.llnl.gov/sites/default/files/inline-files/sundials-2.6.1.tar.gz) and unzip the folder. Relocate the sundials-2.6.2 folder inside the COBRAPRO folder.

3. Download the latest version of [CasADi](https://web.casadi.org/get/) corresponding to your operating system. Unzip and move your CasADi folder inside the COBRAPRO folder. Your COBRAPRO folder should now contain the sundials-2.6.2 and CasADi folders.

4. Before we can install SUNDIALS, the following software are required to compile the mex files that will interface with the SUNDIALS IDA solver:
   - __Mac users__: Download [Xcode](https://developer.apple.com/xcode/) application (can be downloaded from Apple’s App Store). Once Xcode[^1] is installed, proceed to accept the license agreement. This can be done by opening the Xcode application, which will launch a license agreement window and click the “Agree” icon, or type
     ```
     sudo xcodebuild -license accept
     ```
     in Terminal. If the license is not accepted, MATLAB may give an error such as “Xcode is installed, but its license has not been accepted”.
   - __Window users__: Download [MinGW](https://www.mathworks.com/matlabcentral/fileexchange/52848-matlab-support-for-mingw-w64-c-c-fortran-compiler)
     
5. Now you are ready to run `install_COBRAPRO.m`[^2] located inside the COBRAPRO folder. Run `install_COBRAPRO.m` in MATLAB and respond to the prompts displayed in the Command Window in the following manner:

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
Enter return to cancel the installation.
Installation directory:
```
&rarr; Just hit enter
```
Type the name of your CasADi folder (case-sensitive):
```
&rarr; Type the name of the CasADi folder exactly as it appears and hit enter

6. Successful installation will output to the Command Window:
```
COBRAPRO installed successfully! To get started, try running scripts in the Examples folder.
```

[^1]: Note that Xcode requires ~3.4 GB of storage space. 
[^2]: This function will install SUNDIALS by calling the `sundials-2.6.2/ sundialsTB/install_STB.m` file and automatically adds the required folders to your MATLAB path. 
[^3]: COMSOL Multiphsyics is a commerically available finite element analysis software.

## Testing <a name="toc5"></a> ##

Automated test codes are provided in `test` folder:
- `test_1_casadiCheck.m` checks that CasADi is installed and working properly. Successful run will output to Command Window:
```
test_1 successful: CasADi is working properly!
```
- `test_2_comsolValidation.m` validates COBRAPRO against results generated from COMSOL Multiphysics[^3] as a benchmark. This ensures that COBRAPRO is installed properly and that the SUNDIALS IDA solver is working as expected. Successful validation will output to Command Window:
```
test_2 successful: COBRAPRO is working as expected! Results validated against COMSOL.
```

## Examples <a name="toc6"></a> ##

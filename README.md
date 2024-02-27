# COBRAPRO: Co-simulation Battery Modeling for Accelerated Parameter Optimization

<a name="readme-top"></a>

[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]

COBRAPRO is a MATLAB software for physics-based modeling of lithium-ion batteries (LIB) with an embedded parameter identification routine. We aim to provide the battery modeling community with a versatile toolbox for calibrating battery models, a crucial process to achieve accurate simulation results that can capture the behavior of real-world batteries.

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
* MATLAB 2014b and later
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

6. You are ready to use COBRAPRO!

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

<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->



<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/github_username/repo_name">
    <img src="images/logo.png" alt="Logo" width="80" height="80">
  </a>

<h3 align="center">project_title</h3>

  <p align="center">
    project_description
    <br />
    <a href="https://github.com/github_username/repo_name"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/github_username/repo_name">View Demo</a>
    ·
    <a href="https://github.com/github_username/repo_name/issues">Report Bug</a>
    ·
    <a href="https://github.com/github_username/repo_name/issues">Request Feature</a>
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

[![Product Name Screen Shot][product-screenshot]](https://example.com)

Here's a blank template to get started: To avoid retyping too much info. Do a search and replace with your text editor for the following: `github_username`, `repo_name`, `twitter_handle`, `linkedin_username`, `email_client`, `email`, `project_title`, `project_description`

<p align="right">(<a href="#readme-top">back to top</a>)</p>



### Built With

* [![Next][Next.js]][Next-url]
* [![React][React.js]][React-url]
* [![Vue][Vue.js]][Vue-url]
* [![Angular][Angular.io]][Angular-url]
* [![Svelte][Svelte.dev]][Svelte-url]
* [![Laravel][Laravel.com]][Laravel-url]
* [![Bootstrap][Bootstrap.com]][Bootstrap-url]
* [![JQuery][JQuery.com]][JQuery-url]

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- GETTING STARTED -->
## Getting Started

This is an example of how you may give instructions on setting up your project locally.
To get a local copy up and running follow these simple example steps.

### Prerequisites

This is an example of how to list things you need to use the software and how to install them.
* npm
  ```sh
  npm install npm@latest -g
  ```

### Installation

1. Get a free API Key at [https://example.com](https://example.com)
2. Clone the repo
   ```sh
   git clone https://github.com/github_username/repo_name.git
   ```
3. Install NPM packages
   ```sh
   npm install
   ```
4. Enter your API in `config.js`
   ```js
   const API_KEY = 'ENTER YOUR API';
   ```

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- USAGE EXAMPLES -->
## Usage

Use this space to show useful examples of how a project can be used. Additional screenshots, code examples and demos work well in this space. You may also link to more resources.

_For more examples, please refer to the [Documentation](https://example.com)_

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- ROADMAP -->
## Roadmap

- [ ] Feature 1
- [ ] Feature 2
- [ ] Feature 3
    - [ ] Nested Feature

See the [open issues](https://github.com/github_username/repo_name/issues) for a full list of proposed features (and known issues).

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.

If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".
Don't forget to give the project a star! Thanks again!

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTACT -->
## Contact

Your Name - [@twitter_handle](https://twitter.com/twitter_handle) - email@email_client.com

Project Link: [https://github.com/github_username/repo_name](https://github.com/github_username/repo_name)

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

* []()
* []()
* []()

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/github_username/repo_name.svg?style=for-the-badge
[contributors-url]: https://github.com/github_username/repo_name/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/github_username/repo_name.svg?style=for-the-badge
[forks-url]: https://github.com/github_username/repo_name/network/members
[stars-shield]: https://img.shields.io/github/stars/github_username/repo_name.svg?style=for-the-badge
[stars-url]: https://github.com/github_username/repo_name/stargazers
[issues-shield]: https://img.shields.io/github/issues/github_username/repo_name.svg?style=for-the-badge
[issues-url]: https://github.com/github_username/repo_name/issues
[license-shield]: https://img.shields.io/github/license/github_username/repo_name.svg?style=for-the-badge
[license-url]: https://github.com/COBRAPROsimulator/COBRAPRO/blob/main/LICENSE
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/linkedin_username
[product-screenshot]: images/screenshot.png
[Next.js]: https://img.shields.io/badge/next.js-000000?style=for-the-badge&logo=nextdotjs&logoColor=white
[Next-url]: https://nextjs.org/
[React.js]: https://img.shields.io/badge/React-20232A?style=for-the-badge&logo=react&logoColor=61DAFB
[React-url]: https://reactjs.org/
[Vue.js]: https://img.shields.io/badge/Vue.js-35495E?style=for-the-badge&logo=vuedotjs&logoColor=4FC08D
[Vue-url]: https://vuejs.org/
[Angular.io]: https://img.shields.io/badge/Angular-DD0031?style=for-the-badge&logo=angular&logoColor=white
[Angular-url]: https://angular.io/
[Svelte.dev]: https://img.shields.io/badge/Svelte-4A4A55?style=for-the-badge&logo=svelte&logoColor=FF3E00
[Svelte-url]: https://svelte.dev/
[Laravel.com]: https://img.shields.io/badge/Laravel-FF2D20?style=for-the-badge&logo=laravel&logoColor=white
[Laravel-url]: https://laravel.com
[Bootstrap.com]: https://img.shields.io/badge/Bootstrap-563D7C?style=for-the-badge&logo=bootstrap&logoColor=white
[Bootstrap-url]: https://getbootstrap.com
[JQuery.com]: https://img.shields.io/badge/jQuery-0769AD?style=for-the-badge&logo=jquery&logoColor=white
[JQuery-url]: https://jquery.com 

API Documentation
=================

DAE Functions
-------------
.. py:function:: alg_res(x, param)

   Calculates the algebriac residual equations given as g(t,x) = 0, where x includes algebraic and differential variables. The algebraic equations consists of equations for the algebraic variables phis_p, phis_n, phie, jp, and jn.

   :param x: A vector containing the values of algebraic and differential variables.
   :type x: casadi.SX or double
   :param param: A structure containing the model parameters.
   :type param: struct

   :returns: **x_tot** (*casadi.SX or double*) -- A vector where each element corresponds to the solution of the algebraic residual equations g(t,x)=0.

.. py:function:: defineVariableIndex(param)

   Defines the indices and calculates the length for the differential and algebraic variables and stores the information in param structure.

   :param param: A structure containing the model parameters.
   :type param: struct

   :returns: **param** (*struct*) -- A structure containing the model parameters, now inlcuding the variable lengths and indices.

.. py:function:: diff_res(t, x, xp, ida_user_data)

   This function calculates the differential-algebraic equations (DAEs) in implicit form, depending on the method chosen for determining consistent initial conditions:

   * If using the SUNDIALS IDACalcIC method (`param.init_Method = 'IDACalcIC'`), it computes the DAE system F(t, x, x_p) = 0, which includes both the ordinary differential equation (ODE) and algebraic equation (AE) residuals.

   * If using the single-step approach [1] (`param.init_Method = 'SS'`), it calculates the implicit ODE system M(t, x, x_p) = 0 [2], consisting of the perturbed AE and original ODEs with the switch function applied.

   **References**:

   [1] M. T. Lawder, V. Ramadesigan, B. Suthar, and V. R. Subramanian, “Extending explicit and linearly implicit ODE solvers for index-1 DAEs,” Computers & Chemical Engineering, vol. 82, pp. 283–292, Nov. 2015, doi: 10.1016/j.compchemeng.2015.07.002.

   [2] S. Ha and S. Onori, “COBRAPRO: An Open-Source Software for the Doyle-Fuller-Newman Model with Co-Simulation Parameter Optimization Framework,” J. Electrochem. Soc., vol. 171, no. 9, p. 090522, Sep. 2024, doi: 10.1149/1945-7111/ad7292.

   :param t: A structure containing the model parameters.
   :type t: casadi.SX or double
   :param x: A vector containing the values of algebraic and differential variables.
   :type x: casadi.SX or double
   :param xp: A vector containing the values of derivative of the algebraic and differential variables.
   :type xp: casadi.SX or double
   :param ida_user_data: Structure containing additional functions or parameters for the IDA solver.
   :type ida_user_data: struct

   :returns: * **dx_tot** (*casadi.SX or double*) -- Vector comuting F(t, x, x_p) = 0 for the SUNDIALS IDACalcIC method or M(t, x, x_p) = 0 for the single-step approach.
             * **flag** (*double*) -- Required by IDA solver but not used in the code.
             * **new_data** (*double*) -- Required by IDA solver but not used in the code.

.. py:function:: getCurrentDensity(t, param)

   Takes the variable input current profile defined in param.I and calculates the current density at time t.

   :param t: Time instance [s].
   :type t: double
   :param param: A structure containing the model parameters.
   :type param: struct

   :returns: **I** (*double*) -- Current density [A.m-2] at time t.

.. py:function:: jacobianFn(t, x, xp, ~, cj, data)

   Calculates the Jacobian matrix for the discretized DFN model.

   :param t: A structure containing the model parameters.
   :type t: double
   :param x: A vector containing the values of algebraic and differential variables.
   :type x: double
   :param xp: A vector containing the values of derivative of the algebraic and differential variables.
   :type xp: double
   :param cj: Required by IDA to compute Jacobian matrix.
   :type cj: double
   :param data: Structure containing information required to compute Jacobian matrix.
   :type data: struct

   :returns: * **J** (*double*) -- Numerically calculated Jacobian matrix.
             * **flag** (*double*) -- Required by IDA solver but not used in the code.
             * **new_data** (*double*) -- Required by IDA solver but not used in the code.

FVM_interpolation
-----------------

.. py:function:: De_eff_constant(ce, param, domain)

   Calculates the effective electrolyte diffusivity pertaining to the domain of interest when the electrolyte diffusion coefficient is constant. The domain refers to the positive electrode, separator, or negative electrode. 

   :param ce: Normalized electrolyte concentration vector in domain (unitless). 
   :type ce: casadi.SX or double
   :param param: A structure containing the model parameters.
   :type param: struct
   :param domain: Letter indicating battery domain ('p' for positive electrode, 's' for separator, 'n' for negative electrode)
   :type domain: string

   :returns: **Deff** (*casadi.SX or double*) -- Effective electrolyte diffusivity in the domain of interest.

.. py:function:: Kappa_eff_constant(ce, param, domain)

   Calculates the effective electrolyte conductivity pertaining to the domain of interest when the electrolyte conductivity is constant. The domain refers to the positive electrode, separator, or negative electrode. 

   :param ce: Normalized electrolyte concentration vector in domain (unitless). 
   :type ce: casadi.SX or double
   :param param: A structure containing the model parameters.
   :type param: struct
   :param domain: Letter indicating battery domain ('p' for positive electrode, 's' for separator, 'n' for negative electrode)
   :type domain: string

   :returns: **Kappa_eff** (*casadi.SX or double*) -- Effective electrolyte conductivitiy in the domain of interest.

.. py:function:: cs_surf_hermite_interp(cs, r_index_Nr, cs_outer_ghost, param, domain)

   Estimates the surface particle concentration using 3rd order hermite interpolation [1] when using FVM for the radial discretization of the solid particles.

   **References**:

   [1] L. Xu, J. Cooper, A. Allam, and S. Onori, “Comparative Analysis of Numerical Methods for Lithium-Ion Battery Electrochemical Modeling,” J. Electrochem. Soc., vol. 170, no. 12, p. 120525, Dec. 2023, doi: 10.1149/1945-7111/ad1293.

   :param cs: Concatenated normalized solid concentration vector along the particle radius direction for all particles in domain (unitless). 
   :type cs: casadi.SX or double
   :param r_index_Nr: Index of the outer most control volume (closest to particle surface) for the particle of interest.
   :type r_index_Nr: double
   :param cs_outer_ghost: Normalized solid concentration value at the outer ghost node (Nrp+1). 
   :param param: A structure containing the model parameters.
   :type param: struct
   :param domain: Letter indicating battery domain ('p' for positive electrode, 'n' for negative electrode)
   :type domain: string

   :returns: **cs_surf** (*casadi.SX or double*) -- Normalized surface particle concentration.

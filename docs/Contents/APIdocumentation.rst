API Documentation
=================

Functions
---------
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

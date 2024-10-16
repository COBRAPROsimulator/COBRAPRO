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

API Documentation
=================

Functions
---------
.. py:function:: check(molecule, isotope='', linelist='')

   Checks if the parameters passed into the summon() function by the user are valid to use with the ExoMol database.

   :param molecule: Molecule name (e.g. H2O).
   :type molecule: String
   :param isotope: Isotopologue of the molecule (e.g. 1H-16O). The default is ''.
   :type isotope: String, optional
   :param linelist: Line list to be downloaded. The default is ''.
   :type linelist: String, optional

   :rtype: None.

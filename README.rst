Gaia-errors Package
======

**DESCRIPTION:**

Gaia-errors is a Fortran toolkit to apply the Gaia errors to simulated data, 
as specified in the Gaia Science Performance webpapage 



**REQUIREMENTS**

Fortran compile gfortran

**FILES PROVIDED**

- Documentation
   * README.rst

- Tables
   * gfactor-Jun2013.dat: geometrical factors to be applied in the error in parallax as of June 2013
   * TableVr-Oct2014.dat: coefficients of the exponential fit of the error in radial velocity as of October 2014

- Auxialliary routines
   * Compute gaussian dispersion as in Numerical Recipes (gaussdev.f)


Quick Guide
-----------

Given your data in equatorial coordinates (alpha, delta, pi, mu_alpha, mu_delta, Vr), and fixed a population (V,V-I) it returns the errors in astrometry, photometry and atmospherical parameters, and the data affected by errors.

Attribution
-----------

This code is provided by Merce Romero-Gomez - merce.romero at ub.edu, and Josep Manel
Carrasco (UB). The new fit to the parallax error after commissioning have been provided by Rygl, Antoja, DeBruijne et al (2014). The rest of the astrometric errors follow the prescription as in the Gaia Science Performance Webpage (http://www.cosmos.esa.int/web/gaia/science-performance). The atmospheric parameters are implemented using the
values as in Liu & Bailer-Jones (2012). The radial velocity error profile follows the same shape as in the Science Performance webpage and the new fit for the coefficients is
attributed to M.Romero-Gomez, Oct. 2014.

If you have used this code in your research, please let me know and consider acknowledging this package.

License
-------

Copyright (c) 2013-2014 Merce Romero-Gomez

Gaia-errors is open source and free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see `<http://www.gnu.org/licenses/>`_.

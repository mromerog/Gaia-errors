Gaia-errors Package
======

**DESCRIPTION:**

Gaia-errors is a Fortran toolkit to apply the Gaia errors to simulated data, 
as specified in the Gaia Science Performance webpage 



**REQUIREMENTS**

Fortran compile gfortran

**FILES PROVIDED**

- Documentation
   * README.rst

- Tables
   * gfactor-Jun2013.dat: geometrical factors to be applied in the error in parallax as of June 2013
   * TableVr-Jun2015.dat: coefficients of the exponential fit of the error in radial velocity as of October 2014

- Auxialliary routines
   * Compute gaussian dispersion as in Numerical Recipes (gaussdev.f)

- Auxilliary constant files
   * const_math.h: mathematical constants
   * const_ast.h:  astronomical constants
   * const_pop.h:  population constants
   * const_pot.h:  potential constants

Quick Guide
-----------
The main code is main_Gaiaerrors.f. and it calls the routine Gaia-errors to some
input data (provided by the user).
The input file must contain 9 columns with the following:
1-3 cols: x,y,z (in kpc)
4-6 cols: vx,vy,vz (in km/s)
in galactocentric coordinates
col7: Absolute magnitude in V (mag)
col8: (V-I) intrinsic colour (mag)
col9: Absorption in V (mag)
The code performs the change of coordinates from galactocentric to heliocentric and computes the V apparent magnitude and observed (V-I) colour, this is the input necessary for the Gaia-errors subroutine.
The output is transformed back to galactocentric coordinates and it writes input and output in the output file (provided by the user).

To compile the code, just type in the command line:
$ make main_Gaiaerrors
It uses the gfortran compiler.

It also needs, apart from the usual .h files, the cons_pop.h, that includes
a specific Teff, logg and [Fe/H] for a population.

The input data in the Gaia-errors subroutine is:
- month and CAfactor, the length of operational data released in months and the Calibration Astrometric factor to be applied to the astrometric errors, respectively.
- your data in equatorial coordinates (alpha, delta, pi, mu_alpha, mu_delta, Vr)
- a fixed population (V,V-I) 
- a flag regarding the type of errors to be applied (mean end-of-mission or applying geometrical factors and the number of transits)

The output is the errors in astrometry, photometry and atmospherical parameters, and the data affected by errors.

Attribution
-----------

This code is provided by Merce Romero-Gomez - merce.romero at ub.edu, Josep Manel
Carrasco and Roger Mor (UB). The new fit to the parallax error after commissioning have been provided by Rygl, Antoja, DeBruijne et al (2014). The rest of the astrometric errors follow the prescription as in the Gaia Science Performance Webpage (http://www.cosmos.esa.int/web/gaia/science-performance). The atmospheric parameters are implemented using the
values as in Liu & Bailer-Jones (2012). The radial velocity error profile follows the same shape as in the Science Performance webpage.

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

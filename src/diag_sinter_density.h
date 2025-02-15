/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#ifdef DIAG_CLASS
DiagStyle(sinter_density,DiagSinterDensity)

#else

#ifndef SPK_DIAG_SINTER_DENSITY_H
#define SPK_DIAG_SINTER_DENSITY_H

#include "stdio.h"
#include "diag.h"

namespace SPPARKS_NS {

class DiagSinterDensity : public Diag {

 public:
  DiagSinterDensity(class SPPARKS *, int, char **);
  ~DiagSinterDensity() {}
  void init();
  void compute();
  void stats(char *);
  void stats_header(char *);
  
 protected:
  void initialize_parameters_density_calculation(); 

 private:
  class AppSinter *appsinter;
  int nlocal;
  double density, window;
  
  int xstart_density, xend_density, nx;
  int ystart_density, yend_density, ny;
  int zstart_density, zend_density, nz;
};

}

#endif

#endif

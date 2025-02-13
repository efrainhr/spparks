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
DiagStyle(sinter_pore_curvature,DiagSinterPoreCurvature)

#else

#ifndef SPK_DIAG_SINTER_PORE_CURVATURE_H
#define SPK_DIAG_SINTER_PORE_CURVATURE_H

#include "stdio.h"
#include "diag.h"

namespace SPPARKS_NS {

class DiagSinterPoreCurvature : public Diag {
 public:
  DiagSinterPoreCurvature(class SPPARKS *, int, char **);
  ~DiagSinterPoreCurvature() {}
  void init();
  void compute();
  void stats(char *);
  void stats_header(char *);

 private:
  class AppSinter *appsinter;
  double poreCurvature;
  double tripleJunction;	
  double window;
  int **edgeneigh, *faceneigh, surfaceArea;
  void initialize_parameters_calculation();	

	int dimension;

  int xstart_, xend_, nx;
  int ystart_, yend_, ny;
  int zstart_, zend_, nz;
};

}

#endif

#endif

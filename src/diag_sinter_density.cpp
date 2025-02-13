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

#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "output.h"
#include "memory.h"
#include "app.h"
#include "error.h"
#include "timer.h"
#include "diag_sinter_density.h"
#include "app_lattice.h"
#include "app_sinter.h"
#include "comm_lattice.h"
#include "domain.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

DiagSinterDensity::DiagSinterDensity(SPPARKS *spk, int narg, char **arg) : 
  Diag(spk,narg,arg)
{
  if (app->appclass != App::LATTICE)
    error->all(FLERR,"Diag style incompatible with app style");

  if (narg>2)
    error->all(FLERR,"Too many arguments for diag style");

  if (narg==2) {
    window = atof(arg[1]);
		if (window>=0.5)
      error->all(FLERR,"Window size is too large (set to w<0.5)");
	}
  else
    window = 0.333333;
}

/* ---------------------------------------------------------------------- */

void DiagSinterDensity::init()
{
  appsinter = (AppSinter *) app;
  nlocal = appsinter->nlocal;
  density = 0.0; 

  initialize_parameters_density_calculation();
}

/* ---------------------------------------------------------------------- */

void DiagSinterDensity::compute()
{
	double grain_sites = 0;
	double total_sites = 0;
	int xgrid, ygrid, zgrid;
	
	int *spin = appsinter->spin;
	tagint *id = appsinter->id;
	const int VACANT ( AppSinter::VACANT );
	
	for ( int i = 0; i < nlocal; i++ ) {
		appsinter->global_to_grid( id[i], xgrid, ygrid, zgrid );
		if (	((xgrid > xstart_density && xgrid < xend_density) || nx==1) &&
				((ygrid > ystart_density && ygrid < yend_density) || ny==1) &&
				((zgrid > zstart_density && zgrid < zend_density) || nz==1)) {
			total_sites++;
			if ( spin[i] > VACANT ) grain_sites++;
		}
	}

	double local_density_info[2];
  local_density_info[0] = grain_sites;
  local_density_info[1] = total_sites;

	double global_density_info[2] = {0,0};
  MPI_Allreduce(local_density_info, global_density_info, 2, MPI_DOUBLE, MPI_SUM, world);

  density = global_density_info[0] / global_density_info[1];
}

/* ---------------------------------------------------------------------- */

void DiagSinterDensity::stats(char *strtmp)
{
  sprintf(strtmp," %10.6lf",density);
}

/* ---------------------------------------------------------------------- */

void DiagSinterDensity::stats_header(char *strtmp)
{
  sprintf(strtmp," %10s","Density");
}

/* ---------------------------------------------------------------------- */

void DiagSinterDensity::initialize_parameters_density_calculation()
{
  // num sites along each axis in lattice
  nx = domain->nx;
  ny = domain->ny;
  nz = domain->nz;
	 
  xstart_density = static_cast<int>(nx*window);
  xend_density   = static_cast<int>(nx*(1-window));
  ystart_density = static_cast<int>(ny*window);
  yend_density   = static_cast<int>(ny*(1-window));
  zstart_density = static_cast<int>(nz*window);
  zend_density   = static_cast<int>(nz*(1-window));
}

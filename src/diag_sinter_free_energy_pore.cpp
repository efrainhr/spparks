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
#include "diag_sinter_free_energy_pore.h"
#include "app_lattice.h"
#include "app_sinter.h"
#include "comm_lattice.h"
#include "domain.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

DiagSinterFreeEnergyPore::DiagSinterFreeEnergyPore(SPPARKS *spk, int narg, char **arg) : 
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

void DiagSinterFreeEnergyPore::init()
{
  appsinter = (AppSinter *) app;
  nlocal = appsinter->nlocal;
  interfacialFE = 0.0;

  initialize_parameters_calculation();
}

/* ---------------------------------------------------------------------- */

void DiagSinterFreeEnergyPore::compute()
{
	const int VACANT ( AppSinter::VACANT );
	double total_sites = 0;
	int xgrid, ygrid, zgrid;
	
	int *spin = appsinter->spin;
	tagint *id = appsinter->id;
	int *numneigh = appsinter->numneigh;
	int **neighbor = appsinter->neighbor;
	
	double interfacialFEtmp = 0.0;
	for (int i = 0; i < nlocal; i++) {
		appsinter->global_to_grid( id[i], xgrid, ygrid, zgrid );
		if ( ((xgrid > xstart_ && xgrid < xend_) || nx==1) &&
			   ((ygrid > ystart_ && ygrid < yend_) || ny==1) &&
			   ((zgrid > zstart_ && zgrid < zend_) || nz==1) ) {
		  
			total_sites++;
		  
			int ispin = spin[i];

      // If I am a grain site add the number of neighbors that are pore sites
			if ( ispin > VACANT ) {
			  double surface = 0;
				for (int j = 0; j < numneigh[i]; j++)
					if (spin[neighbor[i][j]] == VACANT) surface++;
			  interfacialFEtmp += surface;
			}
		}
	}
	
	vector<double> local_info( 2 );
	local_info[0] = interfacialFEtmp;
	local_info[1] = total_sites;
	
	vector<double> info_all( 2, 0 );
	MPI_Allreduce(&local_info[0], &info_all[0], 3, MPI_DOUBLE, MPI_SUM, world);
	
	interfacialFE = info_all[0] / info_all[1];	
}

/* ---------------------------------------------------------------------- */

void DiagSinterFreeEnergyPore::stats(char *strtmp)
{
  sprintf(strtmp," %10.6lf ",interfacialFE);
}

/* ---------------------------------------------------------------------- */

void DiagSinterFreeEnergyPore::stats_header(char *strtmp)
{
  sprintf(strtmp," %10s ","FE_psite");
}

/* ---------------------------------------------------------------------- */
void DiagSinterFreeEnergyPore::initialize_parameters_calculation()
{
  // num sites along each axis in lattice
  nx = domain->nx;
  ny = domain->ny;
  nz = domain->nz;
  
  xstart_ = static_cast<int>(nx*window);
  xend_   = static_cast<int>(nx*(1-window));
  ystart_ = static_cast<int>(ny*window);
  yend_   = static_cast<int>(ny*(1-window));
  zstart_ = static_cast<int>(nz*window);
  zend_   = static_cast<int>(nz*(1-window));
}

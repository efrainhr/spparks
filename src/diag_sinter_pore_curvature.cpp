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
#include "domain.h"
#include "stdlib.h"
#include "string.h"
#include "output.h"
#include "lattice.h"
#include "memory.h"
#include "app.h"
#include "math.h"
#include "error.h"
#include "timer.h"
#include "diag_sinter_pore_curvature.h"
#include "app_lattice.h"
#include "app_sinter.h"
#include "comm_lattice.h"

using namespace SPPARKS_NS;

// same as in create_sites.cpp and diag_cluster.cpp and lattice.cpp
enum{NONE,LINE_2N,SQ_4N,SQ_8N,TRI,SC_6N,SC_26N,FCC,BCC,DIAMOND,
  FCC_OCTA_TETRA,RANDOM_1D,RANDOM_2D,RANDOM_3D};

/* ---------------------------------------------------------------------- */

DiagSinterPoreCurvature::DiagSinterPoreCurvature(SPPARKS *spk, int narg, char **arg) : 
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

void DiagSinterPoreCurvature::init()
{
  appsinter = (AppSinter *) app;
  poreCurvature = 0.0;
  tripleJunction = 0.0;	

	dimension = domain->dimension;
	int style = domain->lattice->style;

	if (dimension==3 && style==SC_26N) {
    // 12 edges, 3 neighbors per edge
    memory->create(edgeneigh,12,3,"diagcurvature:edgeneigh");

    edgeneigh[0][0] = 21; //(x+1,  y,  z)
    edgeneigh[0][1] = 15; //(  x,y+1,  z)
    edgeneigh[0][2] = 24; //(x+1,y+1,  z)

	  edgeneigh[1][0] = 21; //(x+1,  y,  z)
    edgeneigh[1][1] = 13; //(  x,  y,z+1)
    edgeneigh[1][2] = 22; //(x+1,  y,z+1)

	  edgeneigh[2][0] = 15; //(  x,y+1,  z)
    edgeneigh[2][1] = 13; //(  x,  y,z+1)
    edgeneigh[2][2] = 16; //(  x,y+1,z+1)

	  edgeneigh[3][0] = 15; //(  x,y+1,  z)
    edgeneigh[3][1] = 4;  //(x-1,  y,  z)
    edgeneigh[3][2] = 7;  //(x-1,y+1,  z)

	  edgeneigh[4][0] = 13; //(  x,  y,z+1)
    edgeneigh[4][1] = 4;  //(x-1,  y,  z)
    edgeneigh[4][2] = 5;  //(x-1,  y,z+1)

	  edgeneigh[5][0] = 4;  //(x-1,  y,  z)
    edgeneigh[5][1] = 10; //(  x,y-1,  z)
    edgeneigh[5][2] = 1;  //(x-1,y-1,  z)

	  edgeneigh[6][0] = 4;  //(x-1,  y,  z)
    edgeneigh[6][1] = 12; //(  x,  y,z-1)
    edgeneigh[6][2] = 3;  //(x-1,  y,z-1)

	  edgeneigh[7][0] = 21; //(x+1,  y,  z)
    edgeneigh[7][1] = 10; //(  x,y-1,  z)
    edgeneigh[7][2] = 18; //(x+1,y-1,  z)

	  edgeneigh[8][0] = 13; //(  x,  y,z+1)
    edgeneigh[8][1] = 10; //(  x,y-1,  z)
    edgeneigh[8][2] = 11; //(  x,y-1,z+1)

	  edgeneigh[9][0] = 10; //(  x,y-1,  z)
    edgeneigh[9][1] = 12; //(  x,  y,z-1)
    edgeneigh[9][2] = 9;  //(  x,y-1,z-1)

	  edgeneigh[10][0] = 21; //(x+1,  y,  z)
    edgeneigh[10][1] = 12; //(  x,  y,z-1)
    edgeneigh[10][2] = 20; //(x+1,  y,z-1)

	  edgeneigh[11][0] = 15; //(  x,y+1,  z)
    edgeneigh[11][1] = 12; //(  x,  y,z-1)
    edgeneigh[11][2] = 14; //(  x,y+1,z-1)

    memory->create(faceneigh,6,"diagcurvature:faceneigh");
		faceneigh[0] = 21; //x+1
		faceneigh[1] =  4; //x-1
		faceneigh[2] = 15; //y+1
		faceneigh[3] = 10; //y-1
		faceneigh[4] = 13; //z+1
		faceneigh[5] = 12; //z-1
  }
	else if (dimension==2 && style==SQ_8N) {
    // 4 edges/vertices, 3 neighbors per edge/vertex
    memory->create(edgeneigh,4,3,"diagcurvature:edgeneigh");

    edgeneigh[0][0] = 6; //(x+1,  y)
    edgeneigh[0][1] = 4; //(  x,y+1)
    edgeneigh[0][2] = 7; //(x+1,y+1)

    edgeneigh[1][0] = 1; //(x-1,  y)
    edgeneigh[1][1] = 4; //(  x,y+1)
    edgeneigh[1][2] = 2; //(x-1,y+1)

    edgeneigh[2][0] = 6; //(x+1,  y)
    edgeneigh[2][1] = 3; //(  x,y-1)
    edgeneigh[2][2] = 5; //(x+1,y-1)

	  edgeneigh[3][0] = 1; //(x-1,  y)
    edgeneigh[3][1] = 3; //(  x,y-1)
    edgeneigh[3][2] = 0; //(x-1,y-1)

    memory->create(faceneigh,4,"diagcurvature:faceneigh");
    faceneigh[0] = 6; //x+1
    faceneigh[1] = 1; //x-1
    faceneigh[2] = 3; //y+1
    faceneigh[3] = 4; //y-1
  }
  else
    error->all(FLERR,"Lattice style not compatible with diagnostic");

  initialize_parameters_calculation();
}

/* ---------------------------------------------------------------------- */

void DiagSinterPoreCurvature::compute()
{	
	const int VACANT ( AppSinter::VACANT );
	double pore_sites = 0, saTmp = 0;
	int xgrid, ygrid, zgrid;

	appsinter->comm->all();

	int nlocal = appsinter->nlocal;

	int *spin = appsinter->spin;
	tagint *id = appsinter->id;
	int *numneigh = appsinter->numneigh;
	int **neighbor = appsinter->neighbor;

	double innies = 0, outties = 0, exceptions = 0;
	double tripleJunction_outties = 0, tripleJunction_sameface = 0;

  int nedges = dimension*pow(2,dimension-1);

	for (int i = 0; i < nlocal; i++) {
		appsinter->global_to_grid( id[i], xgrid, ygrid, zgrid );
    if ( ((xgrid > xstart_ && xgrid < xend_) || nx==1) &&
			   ((ygrid > ystart_ && ygrid < yend_) || ny==1) &&
			   ((zgrid > zstart_ && zgrid < zend_) || nz==1) ) {

			int ispin = spin[i], nbor;
			if ( ispin == VACANT ) { // If I am a pore site
				pore_sites++;

				for (int n=0; n<2*dimension; n++)
          if (spin[neighbor[i][faceneigh[n]]] > VACANT) saTmp++;

				for(nbor = 0; nbor < numneigh[i]; nbor++) { // check that it is pore surface
					int neigh = neighbor[i][nbor];
					int value = spin[neigh];
					if ( value > VACANT ) {	// At least one of the neighbors is grain site
						break;
					}
				}
				if ( nbor == numneigh[i] ) { // Previous cycle complete -> no grain site found in the neighborhood
					continue;
				}
				for ( int edge = 0; edge < nedges; edge++ ) {
					int likeENeighs = 0;
					for ( int enbor = 0; enbor < 3; enbor++ ) {
						int eneigh = neighbor[i][edgeneigh[edge][enbor]];
						if ( spin[eneigh] == VACANT ) { // Another pore site
							likeENeighs++;
						}
					}
					if (likeENeighs == 0) {
						outties++;
						// Look for triple junction 
						if ( (spin[neighbor[i][edgeneigh[edge][0]]] != spin[neighbor[i][edgeneigh[edge][1]]])
								|| (spin[neighbor[i][edgeneigh[edge][1]]] != spin[neighbor[i][edgeneigh[edge][2]]]) ) {
							tripleJunction_outties++;
						}
					}
			   	else if(likeENeighs == 2) innies++;
			   	else if(likeENeighs == 1){
				   	if(spin[neighbor[i][edgeneigh[edge][2]]] == VACANT)
							exceptions++;
						else { // Look for triple junction 
							if ( (spin[neighbor[i][edgeneigh[edge][0]]] != spin[neighbor[i][edgeneigh[edge][1]]])
								&& (spin[neighbor[i][edgeneigh[edge][1]]] != spin[neighbor[i][edgeneigh[edge][2]]]) ) {
								tripleJunction_sameface++;
							}
						}
					}
				}
			}
		}
	}
	// exceptions are not considered part of the curvature calculation
	double poreCurvaturetmp = outties - tripleJunction_outties - innies / 3;
	double total_tripleJunction = tripleJunction_outties + tripleJunction_sameface / 2;

	vector<double> local_info( 4 );
	local_info[0] = poreCurvaturetmp;
	local_info[1] = pore_sites;
	local_info[2] = total_tripleJunction;
	local_info[3] = saTmp;

	vector<double> info_all( 4, 0 );
	
	MPI_Allreduce(&local_info[0], &info_all[0], 4, MPI_DOUBLE, MPI_SUM, world);
	
	poreCurvature = info_all[0] / info_all[1];	
	tripleJunction = info_all[2] / info_all[1];
	surfaceArea = static_cast<int>(info_all[3]);
}

/* ---------------------------------------------------------------------- */

void DiagSinterPoreCurvature::stats(char *strtmp)
{
  sprintf(strtmp," %10.6lf %10.6lf %10d ", poreCurvature, tripleJunction, surfaceArea);
}

/* ---------------------------------------------------------------------- */

void DiagSinterPoreCurvature::stats_header(char *strtmp)
{
  sprintf(strtmp," %10s %10s %10s ","PCurvature", "TripleJunc", "SurfaceArea");
}

void DiagSinterPoreCurvature::initialize_parameters_calculation()
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

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

#ifndef MPI_STUBS
#define MPI_STUBS

/* Dummy defs for MPI stubs */

#define MPI_COMM_WORLD 0

#define MPI_SUCCESS 0

#define MPI_INT 1
#define MPI_FLOAT 2
#define MPI_DOUBLE 3
#define MPI_CHAR 4
#define MPI_BYTE 5
#define MPI_LONG_LONG 6
#define MPI_DOUBLE_INT 7

#define MPI_SUM 1
#define MPI_MAX 2
#define MPI_MIN 3
#define MPI_MAXLOC 4
#define MPI_MINLOC 5
#define MPI_LOR 6

#define MPI_ANY_SOURCE -1

#define MPI_Comm int
#define MPI_Request int
#define MPI_Datatype int
#define MPI_Op int

#define MPI_MAX_PROCESSOR_NAME 128

#ifdef __cplusplus
extern "C" {
#endif

/* MPI data structs */

struct MPI_Status {
  int MPI_SOURCE;
};

/* Function prototypes for MPI stubs */

int MPI_Init(int *argc, char ***argv);
int MPI_Initialized(int *flag);
void MPI_Get_processor_name(char *name, int *resultlen);

int MPI_Comm_rank(MPI_Comm comm, int *me);
int MPI_Comm_size(MPI_Comm comm, int *nprocs);
int MPI_Abort(MPI_Comm comm, int errorcode);
int MPI_Finalize();
double MPI_Wtime();

int MPI_Type_size(int, int *);

int MPI_Send(void *buf, int count, MPI_Datatype datatype,
             int dest, int tag, MPI_Comm comm);
int MPI_Rsend(void *buf, int count, MPI_Datatype datatype,
              int dest, int tag, MPI_Comm comm);
int MPI_Recv(void *buf, int count, MPI_Datatype datatype,
             int source, int tag, MPI_Comm comm, MPI_Status *status);
int MPI_Irecv(void *buf, int count, MPI_Datatype datatype,
              int source, int tag, MPI_Comm comm, MPI_Request *request);
int MPI_Wait(MPI_Request *request, MPI_Status *status);
int MPI_Waitall(int n, MPI_Request *request, MPI_Status *status);
int MPI_Waitany(int count, MPI_Request *request, int *index, 
                MPI_Status *status);
int MPI_Sendrecv(void *sbuf, int scount, MPI_Datatype sdatatype,
		  int dest, int stag, void *rbuf, int rcount,
		  MPI_Datatype rdatatype, int source, int rtag,
		  MPI_Comm comm, MPI_Status *status);
int MPI_Get_count(MPI_Status *status, MPI_Datatype datatype, int *count);

int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *comm_out);
int MPI_Comm_dup(MPI_Comm comm, MPI_Comm *comm_out);
int MPI_Comm_free(MPI_Comm *comm);

int MPI_Cart_create(MPI_Comm comm_old, int ndims, int *dims, int *periods,
                    int reorder, MPI_Comm *comm_cart);
int MPI_Cart_get(MPI_Comm comm, int maxdims, int *dims, int *periods,
                 int *coords);
int MPI_Cart_shift(MPI_Comm comm, int direction, int displ,
                   int *source, int *dest);
int MPI_Cart_rank(MPI_Comm comm, int *coords, int *rank);

int MPI_Barrier(MPI_Comm comm);
int MPI_Bcast(void *buf, int count, MPI_Datatype datatype,
              int root, MPI_Comm comm);
int MPI_Allreduce(void *sendbuf, void *recvbuf, int count,
                  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int MPI_Reduce(void *sendbuf, void *recvbuf, int count,
		   MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);
int MPI_Scan(void *sendbuf, void *recvbuf, int count,
             MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int MPI_Allgather(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                  void *recvbuf, int recvcount, MPI_Datatype recvtype,
                  MPI_Comm comm);
int MPI_Allgatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                   void *recvbuf, int *recvcounts, int *displs,
                   MPI_Datatype recvtype, MPI_Comm comm);
int MPI_Reduce_scatter(void *sendbuf, void *recvbuf, int *recvcounts,
                       MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int MPI_Gather(void *sendbuf, int sendcount, MPI_Datatype sendtype,
               void *recvbuf, int recvcount, MPI_Datatype recvtype,
               int root, MPI_Comm comm);
int MPI_Gatherv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
		    void *recvbuf, int *recvcounts, int *displs,
		    MPI_Datatype recvtype, int root, MPI_Comm comm);
int MPI_Alltoallv(const void *sendbuf, const int *sendcounts,
                  const int *sdispls, MPI_Datatype sendtype,
                  void *recvbuf, const int *recvcounts, const int *rdispls,
                  MPI_Datatype recvtype, MPI_Comm comm);

#ifdef __cplusplus
}
#endif

#endif

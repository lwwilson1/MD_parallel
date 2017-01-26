#ifndef H_ROUTINES_H
#define H_ROUTINES_H

#include <mpi.h>
#include "structs.h"

/* Global variable to track particle exchanges */
extern int exchange;

void direct_forces(double *pS,  double *qS,  double *fS,  double *denergy,  int numpar,
                   double *pS2, double *qS2, double *fS2, double *denergy2, int numpar2,
                   int pot_type, double kappa);

void direct_forces_within(struct par *parlist, struct foreng *forlist,
                          int numpar, int pot_type, double kappa);

void perform_ORB(struct par **parlist, int *numparloc, char binrank[20], int numbis,
                 MPI_Comm commarr[numbis], MPI_Datatype MPI_PARTYPE);

void contsruct_commarr(int rank, int numbis, char binrank[20], MPI_Comm **commarr);

void readin_parlist(struct par *parlist, int numparloc, int globparloc, char *sampin1);

void construct_MPI_PARTYPE(MPI_Datatype *MPI_PARTYPE);

void construct_MPI_TNODETYPE(MPI_Datatype *MPI_TNODETYPE);

#endif /* H_ROUTINES_H */

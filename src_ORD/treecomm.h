#ifndef H_TREECOMM_H
#define H_TREECOMM_H

#include <mpi.h>
#include "structs.h"

void gather_all_trees(struct par *parlist, struct tnode *p,
                      struct tnode *globtrees, int wsize,
                      MPI_Datatype MPI_TNODETYPE);

void pc_external_trees(struct tnode *globtrees, struct par *parlist,
                       int numpars, double kappa, int rank, int wsize);

void compute_pc_external(struct tnode *p, double *peng, double *force,
                         double kappa);

#endif /* H_TREECOMM_H */

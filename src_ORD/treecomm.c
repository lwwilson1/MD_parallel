#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

#include "globtreevars.h"

#include "array.h"
#include "tools.h"
#include "treecode.h"
#include "treecomm.h"


void gather_all_trees(struct par *parlist, struct tnode *p,
                      struct tnode *globtrees, int wsize,
                      MPI_Datatype MPI_TNODETYPE)
{
    int i, j, k, l, m, index, offset;
    struct tnodesend psend[9];
    struct tnodesend *precv;
    double *mssend, *msrecv;
    
    int msize = (torder+1)*(torder+1)*(torder+1);
    
    mssend = malloc(sizeof(double) * msize * 9);
    msrecv = malloc(sizeof(double) * msize * 9 * wsize);
    precv  = malloc(sizeof(struct tnodesend) * 9 * wsize);
    
    /* Constructing moments of root node */
    if (p->exist_ms == 0) {
        make_3array(p->ms, torder+1, torder+1, torder+1);
        comp_ms(p, parlist);
        p->exist_ms = 1;
    }
    
    
    /* Packaging root node for sending */
    psend[0].numpar = p->numpar;
    psend[0].ibeg = p->ibeg;
    psend[0].iend = p->iend;
    
    psend[0].x_min = p->x_min;
    psend[0].y_min = p->y_min;
    psend[0].z_min = p->z_min;
    
    psend[0].x_max = p->x_max;
    psend[0].y_max = p->y_max;
    psend[0].z_max = p->z_max;
    
    psend[0].x_mid = p->x_mid;
    psend[0].y_mid = p->y_mid;
    psend[0].z_mid = p->z_mid;
    
    psend[0].radius = p->radius;
    psend[0].sqradius = p->sqradius;
    psend[0].aspect = p->aspect;
    
    psend[0].level = p->level;
    psend[0].num_children = p->num_children;
    psend[0].exist_ms = p->exist_ms;
    psend[0].local = 0;
    
    
    /* Packaging root moments for sending */
    for (i = 0; i < torder + 1; i++) {
        for (j = 0; j < torder + 1; j++) {
            for (k = 0; k < torder + 1; k++) {
                index = i*(torder+1)*(torder+1) + j*(torder+1) + k;
                mssend[index] = p->ms[i][j][k];
            }
        }
    }
    
    for (i = 0; i < p->num_children; i++) {
        
        if (p->child[i]->exist_ms == 0) {
            make_3array(p->child[i]->ms, torder+1, torder+1, torder+1);
            comp_ms(p->child[i], parlist);
            p->child[i]->exist_ms = 1;
        }
        
        psend[i+1].numpar = p->child[i]->numpar;
        psend[i+1].ibeg = p->child[i]->ibeg;
        psend[i+1].iend = p->child[i]->iend;
        
        psend[i+1].x_min = p->child[i]->x_min;
        psend[i+1].y_min = p->child[i]->y_min;
        psend[i+1].z_min = p->child[i]->z_min;
        
        psend[i+1].x_max = p->child[i]->x_max;
        psend[i+1].y_max = p->child[i]->y_max;
        psend[i+1].z_max = p->child[i]->z_max;
        
        psend[i+1].x_mid = p->child[i]->x_mid;
        psend[i+1].y_mid = p->child[i]->y_mid;
        psend[i+1].z_mid = p->child[i]->z_mid;
        
        psend[i+1].radius = p->child[i]->radius;
        psend[i+1].sqradius = p->child[i]->sqradius;
        psend[i+1].aspect = p->child[i]->aspect;
        
        psend[i+1].level = p->child[i]->level;
        psend[i+1].num_children = p->child[i]->num_children;
        psend[i+1].exist_ms = p->child[i]->exist_ms;
        psend[i+1].local = 0;

    
        for (l = 0; l < torder + 1; l++) {
            for (j = 0; j < torder + 1; j++) {
                for (k = 0; k < torder + 1; k++) {
                    index = l*(torder+1)*(torder+1) + j*(torder+1) + k;
                    mssend[msize*(i+1) + index] = p->child[i]->ms[l][j][k];
                }
            }
        }
    }
    
    
    
    /* Gather all the moments and nodes from other procs */
    MPI_Allgather(psend, 9, MPI_TNODETYPE,
                  precv, 9, MPI_TNODETYPE, MPI_COMM_WORLD);
    MPI_Allgather(mssend, msize*9, MPI_DOUBLE,
                  msrecv, msize*9, MPI_DOUBLE, MPI_COMM_WORLD);
    

    /* Set root nodes locally from global moments and nodes */
    for (i = 0; i < wsize; i++) {
        
        (globtrees[i]).numpar = (precv[i*9]).numpar;
        (globtrees[i]).ibeg = (precv[i*9]).ibeg;
        (globtrees[i]).iend = (precv[i*9]).iend;
        
        (globtrees[i]).x_min = (precv[i*9]).x_min;
        (globtrees[i]).y_min = (precv[i*9]).y_min;
        (globtrees[i]).z_min = (precv[i*9]).z_min;
        
        (globtrees[i]).x_max = (precv[i*9]).x_max;
        (globtrees[i]).y_max = (precv[i*9]).y_max;
        (globtrees[i]).z_max = (precv[i*9]).z_max;
        
        (globtrees[i]).x_mid = (precv[i*9]).x_mid;
        (globtrees[i]).y_mid = (precv[i*9]).y_mid;
        (globtrees[i]).z_mid = (precv[i*9]).z_mid;
        
        (globtrees[i]).radius = (precv[i*9]).radius;
        (globtrees[i]).sqradius = (precv[i*9]).sqradius;
        (globtrees[i]).aspect = (precv[i*9]).aspect;
        
        (globtrees[i]).level = (precv[i*9]).level;
        (globtrees[i]).num_children = (precv[i*9]).num_children;
        (globtrees[i]).exist_ms = (precv[i*9]).exist_ms;
        (globtrees[i]).local = (precv[i*9]).local;
        
        for (l = 0; l < torder + 1; l++) {
            for (j = 0; j < torder + 1; j++) {
                for (k = 0; k < torder + 1; k++) {
                    index = l*(torder+1)*(torder+1) + j*(torder+1) + k;
                    (globtrees[i]).ms[l][j][k] = msrecv[i*9*msize + index];
                }
            }
        }
        
        for (m = 0; m < (globtrees[i]).num_children; m++) {
            
            (globtrees[i].child[m])->numpar = (precv[i*9+m+1]).numpar;
            (globtrees[i].child[m])->ibeg = (precv[i*9+m+1]).ibeg;
            (globtrees[i].child[m])->iend = (precv[i*9+m+1]).iend;
            
            (globtrees[i].child[m])->x_min = (precv[i*9+m+1]).x_min;
            (globtrees[i].child[m])->y_min = (precv[i*9+m+1]).y_min;
            (globtrees[i].child[m])->z_min = (precv[i*9+m+1]).z_min;
            
            (globtrees[i].child[m])->x_max = (precv[i*9+m+1]).x_max;
            (globtrees[i].child[m])->y_max = (precv[i*9+m+1]).y_max;
            (globtrees[i].child[m])->z_max = (precv[i*9+m+1]).z_max;
            
            (globtrees[i].child[m])->x_mid = (precv[i*9+m+1]).x_mid;
            (globtrees[i].child[m])->y_mid = (precv[i*9+m+1]).y_mid;
            (globtrees[i].child[m])->z_mid = (precv[i*9+m+1]).z_mid;
            
            (globtrees[i].child[m])->radius = (precv[i*9+m+1]).radius;
            (globtrees[i].child[m])->sqradius = (precv[i*9+m+1]).sqradius;
            (globtrees[i].child[m])->aspect = (precv[i*9+m+1]).aspect;
            
            (globtrees[i].child[m])->level = (precv[i*9+m+1]).level;
            (globtrees[i].child[m])->num_children = (precv[i*9+m+1]).num_children;
            (globtrees[i].child[m])->exist_ms = (precv[i*9+m+1]).exist_ms;
            (globtrees[i].child[m])->local = (precv[i*9+m+1]).local;
            
            offset = i*9*msize + (m+1)*msize;
            
            for (l = 0; l < torder + 1; l++) {
                for (j = 0; j < torder + 1; j++) {
                    for (k = 0; k < torder + 1; k++) {
                        index = l*(torder+1)*(torder+1) + j*(torder+1) + k;
                        (globtrees[i].child[m])->ms[l][j][k] = msrecv[offset + index];
                    }
                }
            }
        }
    }
    
    free(mssend);
    free(msrecv);
    free(precv);
    
    return;
}




void pc_external_trees(struct tnode *globtrees, struct par *parlist,
                       int numpars, double kappa, int rank, int wsize)
{
    /* local variables */
    int i, j, k;
    double peng, force[3];
    
    for (i = 0; i < numpars; i++) {
        
        tarpos[0] = parlist[i].r[0];
        tarpos[1] = parlist[i].r[1];
        tarpos[2] = parlist[i].r[2];
        tarposq = parlist[i].q;
        
        for (j = 0; j < wsize; j++) {
            if (j == rank) continue;
            if (globtrees[j].num_children > 0) {
                for (k = 0; k < globtrees[j].num_children; k++) {
                    peng = 0.0;
                    force[0] = 0.0;
                    force[1] = 0.0;
                    force[2] = 0.0;
            
                    compute_pc_external(globtrees[j].child[k], &peng, force, kappa);
            
                    //UNCOMMENT THESE FOR OUT OF CELL INTERACTIONS
                    parlist[i].peng += peng * parlist[i].q;
                    parlist[i].f[0] += force[0] * parlist[i].q;
                    parlist[i].f[1] += force[1] * parlist[i].q;
                    parlist[i].f[2] += force[2] * parlist[i].q;
                }
            } else {
                peng = 0.0;
                force[0] = 0.0;
                force[1] = 0.0;
                force[2] = 0.0;
            
                compute_pc_external(&(globtrees[j]), &peng, force, kappa);
            
                //UNCOMMENT THESE FOR OUT OF CELL INTERACTIONS
                parlist[i].peng += peng * parlist[i].q;
                parlist[i].f[0] += force[0] * parlist[i].q;
                parlist[i].f[1] += force[1] * parlist[i].q;
                parlist[i].f[2] += force[2] * parlist[i].q;
            }
        }
    }
    
    return;
}




void compute_pc_external(struct tnode *p, double *peng, double *force,
                         double kappa)
{
    /* local variables */
    double tx, ty, tz;
    int i, j, k;
    
    
    /* determine DISTSQ for MAC test */
    tx = tarpos[0] - p->x_mid;
    ty = tarpos[1] - p->y_mid;
    tz = tarpos[2] - p->z_mid;
    
    
    /* initialize local potential energy and force */
    *peng = 0.0;
    force[0] = 0.0;
    force[1] = 0.0;
    force[2] = 0.0;
    
    /* Compute as if MAC is accepted */
        
    for (i = 0; i < torderlim + 3; i++) {
        for (j = 0; j < torderlim + 3; j++) {
            for (k = 0; k < torderlim + 3; k++) {
                b1[i][j][k] = 0.0;
                a1[i][j][k] = 0.0;
            }
        }
    }
        
    comp_tcoeff_yuk(tx, ty, tz, kappa);
    
    for (k = 0; k < torder+1; k++) {
        for (j = 0; j < torder-k+1; j++) {
            for (i = 0; i < torder-k-j+1; i++) {
                *peng += a1[i+2][j+2][k+2] * p->ms[i][j][k];
                force[0] += cf[i] * a1[i+3][j+2][k+2] * p->ms[i][j][k];
                force[1] += cf[j] * a1[i+2][j+3][k+2] * p->ms[i][j][k];
                force[2] += cf[k] * a1[i+2][j+2][k+3] * p->ms[i][j][k];
            }
        }
    }

    return;
}

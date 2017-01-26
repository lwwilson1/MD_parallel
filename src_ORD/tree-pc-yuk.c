/*
 *Procedures for Cluster Particle Treecode
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "globtreevars.h"

#include "array.h"
#include "tools.h"
#include "structs.h"
#include "treecode.h"


/* Definition of variables declared extern in globtreevars.h */
int torder, torderlim;
double *cf, *cf1, *cf2, *cf3;
double ***a1, ***b1;

int minlevel, maxlevel;

int orderoffset;
double tarpos[3];
double thetasq, tarposq;

int *orderarr;



void setup_yuk(struct par *parlist, int numpars, int order,
               double *xyzminmax, double theta)
{
    /*These are local variables. Definitions not necessary right here.*/
    int i;
    double t1;

    /*changing values of our extern variables*/
    torder = order;
    torderlim = torder + 1;
    thetasq = theta*theta;
    
    /*setting up extern variables for tree level tracking */
    minlevel = 50000;
    maxlevel = 0;
    

    /*allocating global Taylor expansion variables*/
    make_vector(cf, torderlim);
    make_vector(cf1, torderlim);
    make_vector(cf2, torderlim);
    make_vector(cf3, torderlim);

    make_3array(a1, torderlim+3, torderlim+3, torderlim+3);
    make_3array(b1, torderlim+3, torderlim+3, torderlim+3);


    /*initializing arrays for Taylor sums and coefficients*/
    for (i = 0; i < torderlim; i++) {
        cf[i] = i + 1.0;
        t1 = 1.0 / (i + 1.0);
        cf1[i] = t1;
        cf2[i] = 1.0 - (0.5 * t1);
        cf3[i] = 1.0 - t1;
    }

    
    /*find bounds of Cartesian box enclosing the particles*/
    minmaxpar(parlist, numpars, xyzminmax);

    
    /*make and init orderarr to keep up with shuffled particles */
    make_vector(orderarr, numpars);
    for (i = 0; i < numpars; i++) orderarr[i] = i+1;

    
    return;
}


/* TREE_COMPFP routine */
void pc_treecode_yuk(struct tnode *p, struct par *parlist,
                     int numpars, double kappa)
{
    /* local variables */
    int i, j;
    double tempx, tempq, peng, force[3];

    for (i = 0; i < numpars; i++) {
        parlist[i].peng = 0.0;
        parlist[i].f[0] = 0.0;
        parlist[i].f[1] = 0.0;
        parlist[i].f[2] = 0.0;
        
        tarpos[0] = parlist[i].r[0];
        tarpos[1] = parlist[i].r[1];
        tarpos[2] = parlist[i].r[2];
        tarposq = parlist[i].q;

        /* Temporarily moving particle i to prevent self-interaction */
        tempx = parlist[i].r[0];
        tempq = parlist[i].q;
        
        parlist[i].r[0] += 1000.0;
        parlist[i].q = 0.0;

        for (j = 0; j < p->num_children; j++) {
            peng = 0.0;
            force[0] = 0.0;
            force[1] = 0.0;
            force[2] = 0.0;
            
            compute_pc_yuk(p->child[j], parlist, &peng, force,
                           numpars, kappa);
            
            parlist[i].peng += peng;
            parlist[i].f[0] += force[0];
            parlist[i].f[1] += force[1];
            parlist[i].f[2] += force[2];
        }
        
        parlist[i].r[0] = tempx;
        parlist[i].q = tempq;
        
        parlist[i].peng *= tempq;
        parlist[i].f[0] *= tempq;
        parlist[i].f[1] *= tempq;
        parlist[i].f[2] *= tempq;
    }

    return;

} /* END of function pc_treecode_yuk */



/* COMPFP_TREE routine */
void compute_pc_yuk(struct tnode *p, struct par *parlist,
                    double *peng, double *force, int arrdim, double kappa)
{
    /* local variables */
    double tx, ty, tz, distsq;
    double penglocal, flocal[3];
    int i, j, k;


    /* determine DISTSQ for MAC test */
    tx = tarpos[0] - p->x_mid;
    ty = tarpos[1] - p->y_mid;
    tz = tarpos[2] - p->z_mid;
    distsq = tx*tx + ty*ty + tz*tz;
    

    /* initialize local potential energy and force */
    *peng = 0.0;
    force[0] = 0.0;
    force[1] = 0.0;
    force[2] = 0.0;
    

    if ((p->sqradius < distsq * thetasq) && (p->sqradius != 0.00)) {

    /*
     * If MAC is accepted and there is more than 1 particle
     * in the box, use the expansion for the approximation.
     */
    
        for (i = 0; i < torderlim + 3; i++) {
            for (j = 0; j < torderlim + 3; j++) {
                for (k = 0; k < torderlim + 3; k++) {
                    b1[i][j][k] = 0.0;
                    a1[i][j][k] = 0.0;
                }
            }
        }

        if (p->exist_ms == 0) {
            make_3array(p->ms, torder+1, torder+1, torder+1);
            comp_ms(p, parlist);
            p->exist_ms = 1;
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

    } else {

    /*
     * If MAC fails check to see if there are children. If not, perform direct
     * calculation. If there are children, call routine recursively for each.
     */
        
        if (p->num_children == 0) {
            comp_direct_yuk(&penglocal, flocal, p->ibeg, p->iend,
                            parlist, kappa);
            
            *peng = penglocal;
            force[0] = flocal[0];
            force[1] = flocal[1];
            force[2] = flocal[2];
        } else {
            for (i = 0; i < p->num_children; i++) {
                compute_pc_yuk(p->child[i], parlist, &penglocal, flocal,
                               arrdim, kappa);
                
                *peng += penglocal;
                force[0] += flocal[0];
                force[1] += flocal[1];
                force[2] += flocal[2];
            }
        }
    }


    return;

} /* END of function compute_cp1 */




void comp_tcoeff_yuk(double dx, double dy, double dz, double kappa)
{
        /* local variables */
        double fac, sqfac, dist, dist2;
        double tdx, tdy, tdz, kappax, kappay, kappaz;
        int i, j, k, i1, i2, j1, j2, k1, k2;


        /* setup variables */
        tdx = 2.0 * dx;
        tdy = 2.0 * dy;
        tdz = 2.0 * dz;

        kappax = kappa * dx;
        kappay = kappa * dy;
        kappaz = kappa * dz;

        dist2 = dx*dx + dy*dy + dz*dz;
        fac = 1.0 / (dist2);
        sqfac = sqrt(fac);
        dist = sqrt(dist2);


        /* 0th coefficient or function val */
        a1[0][0][0] = exp(-kappa * dist);
        b1[0][0][0] = a1[0][0][0] / dist;


        /* set of indices for which two of them are 0 */
        a1[1][0][0] = kappax * b1[0][0][0];
        a1[0][1][0] = kappay * b1[0][0][0];
        a1[0][0][1] = kappaz * b1[0][0][0];

        b1[1][0][0] = fac * dx * (b1[0][0][0] + kappa*a1[0][0][0]);
        b1[0][1][0] = fac * dy * (b1[0][0][0] + kappa*a1[0][0][0]);
        b1[0][0][1] = fac * dz * (b1[0][0][0] + kappa*a1[0][0][0]);

        for (i = 2; i < torderlim+1; i++) {
                i1 = i - 1;
                i2 = i - 2;

                a1[i][0][0] = cf1[i1] * kappa * (dx * b1[i1][0][0] - b1[i2][0][0]);
                a1[0][i][0] = cf1[i1] * kappa * (dy * b1[0][i1][0] - b1[0][i2][0]);
                a1[0][0][i] = cf1[i1] * kappa * (dz * b1[0][0][i2] - b1[0][0][i2]);

                b1[i][0][0] = fac * (tdx * cf2[i1] * b1[i1][0][0]
                                         - cf3[i1] * b1[i2][0][0]
                           + cf1[i1] * kappa * (dx * a1[i1][0][0]
                                                   - a1[i2][0][0]));

                b1[0][i][0] = fac * (tdy * cf2[i1] * b1[0][i1][0]
                                         - cf3[i1] * b1[0][i2][0]
                           + cf1[i1] * kappa * (dy * a1[0][i1][0]
                                                   - a1[0][i2][0]));

                b1[0][0][i] = fac * (tdz * cf2[i1] * b1[0][0][i1]
                                         - cf3[i1] * b1[0][0][i2]
                           + cf1[i1] * kappa * (dz * a1[0][0][i1]
                                                   - a1[0][0][i2]));
        }


        /* set of indices for which one is 0, one is 1, and other is >= 1 */
        a1[1][1][0] = kappax * b1[0][1][0];
        a1[1][0][1] = kappay * b1[0][0][1];
        a1[0][1][1] = kappaz * b1[0][0][1];


        b1[1][1][0] = fac * (dx * b1[0][1][0] + tdy * b1[1][0][0]
                                           + kappax * a1[0][1][0]);
        b1[1][0][1] = fac * (dx * b1[0][0][1] + tdz * b1[1][0][0]
                                           + kappax * a1[0][0][1]);
        b1[0][1][1] = fac * (dy * b1[0][0][1] + tdz * b1[0][1][0]
                                           + kappay * a1[0][0][1]);

        for (i = 2; i < torderlim; i++) {
                i1 = i - 1;
                i2 = i - 2;

                a1[1][0][i] = kappax * b1[0][0][i];
                a1[0][1][i] = kappay * b1[0][0][i];
                a1[0][i][1] = kappaz * b1[0][i][0];
                a1[1][i][0] = kappax * b1[0][i][0];
                a1[i][1][0] = kappay * b1[i][0][0];
                a1[i][0][1] = kappaz * b1[i][0][0];

                b1[1][0][i] = fac * (dx * b1[0][0][i] + tdz * b1[1][0][i1]
                                                            - b1[1][0][i2]
                                                   + kappax * a1[0][0][i]);

                b1[0][1][i] = fac * (dy * b1[0][0][i] + tdz * b1[0][1][i1]
                                        - b1[0][1][i2]
                                                   + kappay * a1[0][0][i]);

                b1[0][i][1] = fac * (dz * b1[0][i][0] + tdy * b1[0][i1][1]
                                        - b1[0][i2][1]
                                                   + kappaz * a1[0][i][0]);

                b1[1][i][0] = fac * (dx * b1[0][i][0] + tdy * b1[1][i1][0]
                                        - b1[1][i2][0]
                                                   + kappax * a1[0][i][0]);

                b1[i][1][0] = fac * (dy * b1[i][0][0] + tdx * b1[i1][1][0]
                                        - b1[i2][1][0]
                                                   + kappay * a1[i][0][0]);

                b1[i][0][1] = fac * (dz * b1[i][0][0] + tdx * b1[i1][0][1]
                                        - b1[i2][0][1]
                                                   + kappaz * a1[i][0][0]);
        }

        /* set of indices for which one is 0, others are >=2 */
        for (i = 2; i < torderlim - 1; i++) {
                i1 = i - 1;
                i2 = i - 2;

                for (j = 2; j < torderlim - i + 1; j++) {
                        j1 = j - 1;
                        j2 = j - 2;

                        a1[i][j][0] = cf1[i1] * kappa * (dx * b1[i1][j][0]
                                                            - b1[i2][j][0]);

                        a1[i][0][j] = cf1[i1] * kappa * (dx * b1[i1][0][j]
                                                            - b1[i2][0][j]);

                        a1[0][i][j] = cf1[i1] * kappa * (dy * b1[0][i1][j]
                                                            - b1[0][i2][j]);
                        
                        b1[i][j][0] = fac * (tdx * cf2[i1] * b1[i1][j][0]
                                                     + tdy * b1[i][j1][0]
                                                 - cf3[i1] * b1[i2][j][0]
                                                           - b1[i][j2][0]
                                   + cf1[i1] * kappa * (dx * a1[i1][j][0]
                                                           - a1[i2][j][0]));

                        b1[i][0][j] = fac * (tdx * cf2[i1] * b1[i1][0][j]
                                                     + tdz * b1[i][0][j1]
                                                 - cf3[i1] * b1[i2][0][j]
                                                           - b1[i][0][j2]
                                   + cf1[i1] * kappa * (dx * a1[i1][0][j]
                                                           - a1[i2][0][j]));

                        b1[0][i][j] = fac * (tdy * cf2[i1] * b1[0][i1][j]
                                                     + tdz * b1[0][i][j1]
                                                 - cf3[i1] * b1[0][i2][j]
                                                           - b1[0][i][j2]
                                   + cf1[i1] * kappa * (dy * a1[0][i1][j]
                                                           - a1[0][i2][j]));
                }
        }

        /* set of indices for which two are 1, other is >= 1 */
        a1[1][1][1] = kappax * b1[0][1][1];
        b1[1][1][1] = fac * (dx * b1[0][1][1] + tdy * b1[1][0][1] 
                                              + tdz * b1[1][1][0]
                                           + kappax * a1[0][1][1]);

        for (i = 2; i < torderlim - 1; i++) {
                i1 = i - 1;
                i2 = i - 2;

                a1[1][1][i] = kappax * b1[0][1][i];
                a1[1][i][1] = kappax * b1[0][i][1];
                a1[i][1][1] = kappay * b1[i][0][1];

                b1[1][1][i] = fac * (dx * b1[0][1][i] + tdy * b1[1][0][i]
                                  + tdz * b1[1][1][i1]      - b1[1][1][i2]
                                                   + kappax * a1[0][1][i]);

                b1[1][i][1] = fac * (dx * b1[0][i][1] + tdy * b1[1][i1][1]
                                  + tdz * b1[1][i][0]       - b1[1][i2][1]
                                                   + kappax * a1[0][i][1]);

                b1[i][1][1] = fac * (dy * b1[i][0][1] + tdx * b1[i1][1][1]
                                  + tdz * b1[i][1][0]       - b1[i2][1][1]
                                                   + kappay * a1[i][0][1]);
        }

        /* set of indices for which one is 1, others are >= 2 */
        for (i = 2; i < torderlim - 2; i++) {
                i1 = i - 1;
                i2 = i - 2;

                for (j = 2; j < torderlim - i + 1; j++) {
                        j1 = j - 1;
                        j2 = j - 2;

                        a1[1][i][j] = kappax * b1[0][i][j];
                        a1[i][1][j] = kappay * b1[i][0][j];
                        a1[i][j][1] = kappaz * b1[i][j][0];

                        b1[1][i][j] = fac * (dx * b1[0][i][j] + tdy * b1[1][i1][j]
                                          + tdz * b1[1][i][j1]      - b1[1][i2][j]
                                                - b1[1][i][j2]
                                                           + kappax * a1[0][i][j]);

                        b1[i][1][j] = fac * (dy * b1[i][0][j] + tdx * b1[i1][1][j]
                                          + tdz * b1[i][1][j1]      - b1[i2][1][j]
                                                - b1[i][1][j2]
                                                           + kappay * a1[i][0][j]);

                        b1[i][j][1] = fac * (dz * b1[i][j][0] + tdx * b1[i1][j][1]
                                          + tdy * b1[i][j1][1]      - b1[i2][j][1]
                                                - b1[i][j2][1]
                                                           + kappaz * a1[i][j][0]);
                }
        }

        /* set of indices for which all are >= 2 */

        //reorder this to make it more friendly to row-major storage
        for (k = 2; k < torderlim - 3; k++) {
                k1 = k - 1;
                k2 = k - 2;
                
                for (j = 2; j < torderlim - k - 1; j++) {
                        j1 = j - 1;
                        j2 = j - 2;

                        for (i = 2; i < torderlim - k - j + 1; i++) {
                                i1 = i - 1;
                                i2 = i - 2;

                                a1[i][j][k] = cf1[i1] * kappa * (dx * b1[i1][j][k]
                                                                    - b1[i2][j][k]);

                                b1[i][j][k] = fac * (tdx * cf2[i1] * b1[i1][j][k]
                                                             + tdy * b1[i][j1][k] 
                                                             + tdz * b1[i][j][k1]
                                                         - cf3[i1] * b1[i2][j][k]
                                                    - b1[i][j2][k] - b1[i][j][k2]
                                           + cf1[i1] * kappa * (dx * a1[i1][j][k]
                                                                   - a1[i2][j][k]));
                        }
                }
        }

        return;

} /* END function comp_tcoeff */




/*
 * comp_direct directly computes the potential on the targets in the current cluster due
 * to the current source, determined by the global variable TARPOS
 */
void comp_direct_yuk(double *peng, double *force, int ibeg, int iend,
                     struct par *parlist, double kappa)
{

    /* local variables */
    int i;
    double tx, ty, tz, dist;
    double temp, qloc;

    *peng = 0.0;
    force[0] = 0.0;
    force[1] = 0.0;
    force[2] = 0.0;
    
    for (i = ibeg - 1; i < iend; i++) {

        tx = parlist[i].r[0] - tarpos[0];
        ty = parlist[i].r[1] - tarpos[1];
        tz = parlist[i].r[2] - tarpos[2];
        dist = sqrt(tx*tx + ty*ty + tz*tz);
        qloc = parlist[i].q;
        
        temp = exp(-kappa * dist) / dist;
        *peng += qloc * temp;
        
        temp *= (kappa * dist + 1.0) / pow(dist, 2);
        force[0] -= qloc * tx * temp;
        force[1] -= qloc * ty * temp;
        force[2] -= qloc * tz * temp;
    }

    return;

} /* END function comp_direct */

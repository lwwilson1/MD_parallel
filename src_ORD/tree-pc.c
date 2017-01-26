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
#include "partition.h"
#include "treecode.h"


void create_tree_n0(struct tnode **p, int ibeg, int iend, struct par *parlist,
                    int maxparnode, double *xyzmm, int level, int arrdim)
{
        /*local variables*/
        double x_mid, y_mid, z_mid, xl, yl, zl, lmax, t1, t2, t3;
        int ind[8][2];
        double xyzmms[6][8];
        int i, j, loclev, numposchild;
        double lxyzmm[6];


        for (i = 0; i < 8; i++)
                for (j = 0; j < 2; j++)
                        ind[i][j] = 0.0;

        for (i = 0; i < 6; i++)
                for (j = 0; j < 8; j++)
                        xyzmms[i][j] = 0.0;

        for (i = 0; i < 6; i++) lxyzmm[i] = 0.0;
                        

        (*p) = malloc(sizeof(struct tnode));


        /* set node fields: number of particles, exist_ms, and xyz bounds */
        (*p)->numpar = iend - ibeg + 1;
        (*p)->exist_ms = 0;


        (*p)->x_min = xyzmm[0];
        (*p)->x_max = xyzmm[1];
        (*p)->y_min = xyzmm[2];
        (*p)->y_max = xyzmm[3];
        (*p)->z_min = xyzmm[4];
        (*p)->z_max = xyzmm[5];



        /*compute aspect ratio*/
        xl = (*p)->x_max - (*p)->x_min;
        yl = (*p)->y_max - (*p)->y_min;
        zl = (*p)->z_max - (*p)->z_min;
        
        lmax = max3(xl, yl, zl);
        t1 = lmax;
        t2 = min3(xl, yl, zl);


        if (t2 != 0.0)
                (*p)->aspect = t1/t2;
        else
                (*p)->aspect = 0.0;


        /*midpoint coordinates, RADIUS and SQRADIUS*/
        (*p)->x_mid = ((*p)->x_max + (*p)->x_min) / 2.0;
        (*p)->y_mid = ((*p)->y_max + (*p)->y_min) / 2.0;
        (*p)->z_mid = ((*p)->z_max + (*p)->z_min) / 2.0;

        t1 = (*p)->x_max - (*p)->x_mid;
        t2 = (*p)->y_max - (*p)->y_mid;
        t3 = (*p)->z_max - (*p)->z_mid;

        (*p)->sqradius = t1*t1 + t2*t2 + t3*t3;
        (*p)->radius = sqrt((*p)->sqradius);

        /*set particle limits, tree level of node, and nullify child pointers*/
        (*p)->ibeg = ibeg;
        (*p)->iend = iend;
        (*p)->level = level;


        if (maxlevel < level)
                maxlevel = level;

        (*p)->num_children = 0;
        for (i = 0; i < 8; i++)
                (*p)->child[i] = NULL;


        if ((*p)->numpar > maxparnode) {

        /*
         * set IND array to 0, and then call PARTITION_8 routine.
         * IND array holds indices of the eight new subregions.
         * Also, setup XYZMMS array in the case that SHRINK = 1.
         */
                xyzmms[0][0] = (*p)->x_min;
                xyzmms[1][0] = (*p)->x_max;
                xyzmms[2][0] = (*p)->y_min;
                xyzmms[3][0] = (*p)->y_max;
                xyzmms[4][0] = (*p)->z_min;
                xyzmms[5][0] = (*p)->z_max;

                ind[0][0] = ibeg;
                ind[0][1] = iend;

                x_mid = (*p)->x_mid;
                y_mid = (*p)->y_mid;
                z_mid = (*p)->z_mid;

                partition_8(parlist, xyzmms, xl, yl, zl, lmax, &numposchild,
                            x_mid, y_mid, z_mid, ind);

                loclev = level + 1;

                for (i = 0; i < numposchild; i++) {
                        if (ind[i][0] <= ind[i][1]) {
                                (*p)->num_children = (*p)->num_children + 1;

                                for (j = 0; j < 6; j++) lxyzmm[j] = xyzmms[j][i];

                                create_tree_n0(
                                        &((*p)->child[(*p)->num_children - 1]),
                                        ind[i][0], ind[i][1], parlist,
                                        maxparnode, lxyzmm, loclev, arrdim);
                        }
                }

        } else {
                if (level < minlevel) minlevel = level;
        }

        return;

} /* end of function create_tree_n0 */





void partition_8(struct par *parlist,
                 double xyzmms[6][8], double xl, double yl, double zl,
                 double lmax, int *numposchild,
                 double x_mid, double y_mid, double z_mid,
                 int ind[8][2])
//IN THE FORTRAN, numposchild is INOUT! I may need to make this a pointer instead
//Note: I'm passing the address from the calls to partition_8
{

        /* local variables */
        int temp_ind, i, j;
        double critlen;

        *numposchild = 1;
        critlen = lmax / sqrt(2.0);


        if (xl >= critlen) {

                partition(parlist, 0, orderarr, ind[0][0], ind[0][1],
                          x_mid, &temp_ind);


                ind[1][0] = temp_ind + 1;
                ind[1][1] = ind[0][1];
                ind[0][1] = temp_ind;

                for (i = 0; i < 6; i++)
                        xyzmms[i][1] = xyzmms[i][0];

                xyzmms[1][0] = x_mid;
                xyzmms[0][1] = x_mid;
                *numposchild = 2 * *numposchild;

        }


        if (yl >= critlen) {

                for (i = 0; i < *numposchild; i++) {
                        partition(parlist, 1, orderarr, ind[i][0], ind[i][1],
                                  y_mid, &temp_ind);
                        
                        ind[*numposchild + i][0] = temp_ind + 1;
                        ind[*numposchild + i][1] = ind[i][1];
                        ind[i][1] = temp_ind;

                        for (j = 0; j < 6; j++)
                                xyzmms[j][*numposchild + i] = xyzmms[j][i];

                        xyzmms[3][i] = y_mid;
                        xyzmms[2][*numposchild + i] = y_mid;
                }

                *numposchild = 2 * *numposchild;

        }


        if (zl >= critlen) {

                for (i = 0; i < *numposchild; i++) {
                        partition(parlist, 2, orderarr, ind[i][0], ind[i][1],
                                  z_mid, &temp_ind);
                        
                        ind[*numposchild + i][0] = temp_ind + 1;
                        ind[*numposchild + i][1] = ind[i][1];
                        ind[i][1] = temp_ind;

                        for (j = 0; j < 6; j++)
                                xyzmms[j][*numposchild + i] = xyzmms[j][i];

                        xyzmms[5][i] = z_mid;
                        xyzmms[4][*numposchild + i] = z_mid;
                }

                *numposchild = 2 * *numposchild;

        }

        return;

} /* END of function partition_8 */






/*
 * comp_ms computes the moments for node p needed in the Taylor approximation
 */
void comp_ms(struct tnode *p, struct par *parlist)
{
    int i, k1, k2, k3;
    double dx, dy, dz, tx, ty, tz, qloc;
    
    for (k1 = 0; k1 < torder + 1; k1++) {
        for (k2 = 0; k2 < torder + 1; k2++) {
            for (k3 = 0; k3 < torder + 1; k3++) {
                p->ms[k1][k2][k3] = 0.0;
            }
        }
    }
    
    for (i = p->ibeg-1; i < p->iend; i++) {
        dx = parlist[i].r[0] - p->x_mid;
        dy = parlist[i].r[1] - p->y_mid;
        dz = parlist[i].r[2] - p->z_mid;
        qloc = parlist[i].q;

        tz = 1.0;
        for (k3 = 0; k3 < torder + 1; k3++) {
            ty = 1.0;
            for (k2 = 0; k2 < torder - k3 + 1; k2++) {
                tx = 1.0;
                for (k1 = 0; k1 < torder - k3 - k2 + 1; k1++) {
                    p->ms[k1][k2][k3] += qloc * tx*ty*tz;
                    tx *= dx;
                }
                ty *= dy;
            }
            tz *= dz;
        }
    }

    return;

} /* END function comp_ms */






/*
 * cleanup deallocates allocated global variables and then calls recursive function
 * remove_node to delete the tree
 */
void cleanup(struct tnode *p)
{
    free_vector(cf);
    free_vector(cf1);
    free_vector(cf2);
    free_vector(cf3);
    free_3array(b1);
    free_3array(a1);

    free_vector(orderarr);
    
    remove_node(p);
    free(p);

    return;

} /* END function cleanup */




void remove_node(struct tnode *p)
{
    /* local variables */
    int i;

    if (p->exist_ms == 1)
            free(p->ms);

    if (p->num_children > 0) {
        for (i = 0; i < p->num_children; i++) {
            remove_node(p->child[i]);
            free(p->child[i]);
        }
    }

    return;

} /* END function remove_node */

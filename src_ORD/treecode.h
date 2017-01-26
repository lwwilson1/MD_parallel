#ifndef H_TREEFUNCTIONS_H
#define H_TREEFUNCTIONS_H

#include "structs.h"

/* declaration of treecode support functions*/

void remove_node(struct tnode *p);

void cleanup(struct tnode *p);

void comp_direct_yuk(double *peng, double *force, int ibeg, int iend,
                     struct par *parlist, double kappa);

void comp_ms(struct tnode *p, struct par *parlist);

void comp_tcoeff_yuk(double dx, double dy, double dz,
                     double kappa);

void compute_pc_yuk(struct tnode *p, struct par *parlist,
                    double *peng, double *force, int arrdim, double kappa);

void pc_treecode_yuk(struct tnode *p, struct par *parlist,
                     int numpars, double kappa);

void partition_8(struct par *parlist,
                 double xyzmms[6][8], double xl, double yl, double zl,
                 double lmax, int *numposchild,
                 double x_mid, double y_mid, double z_mid,
                 int ind[8][2]);

void create_tree_n0(struct tnode **p, int ibeg, int iend,
                    struct par *parlist,
                    int maxparnode, double *xyzmm,
                    int level, int arrdim);

void setup_yuk(struct par *parlist, int numpars, int order,
               double *xyzminmax, double theta);

#endif /* H_TREEFUNCTIONS_H */

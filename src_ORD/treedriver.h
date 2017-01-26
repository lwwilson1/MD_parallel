#ifndef H_TREEDRIVER_H
#define H_TREEDRIVER_H

#include "structs.h"

/* declaration of primary treecode driver */

void treecode(struct par *parlist, int numpars, struct foreng *forlist,
              int order, double theta, int maxparnode, int *orderind,
              double *timetree, int pot_type, double kappa);

#endif /* H_TREEFUNCTIONS_H */

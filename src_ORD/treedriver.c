#include <time.h>
//#include <sys/time.h>
#include <stdio.h>

#include "array.h"
#include "globvars.h"
#include "tnode.h"
#include "tools.h"
#include "tree-pc.h"

#include "treedriver.h"
/* definition of primary treecode driver */

// particle has parlist.r[0], .r[1], .r[2], .q for *xS, *yS, *zS, *qS
// foreng has forlist.f[0], .f[1], .f[2], .pengy for *tpeng, *tforce

void treecode(struct par *parlist, int numpars, struct foreng *forlist,
              int order, double theta, int maxparnode, int *orderind,
              double *timetree, int pot_type, double kappa)
{

        /* local variables */
        struct tnode *troot;
        int i, level;
        double xyzminmax[6];
    

        /* date and time */
        time_t time1, time2;
        double totaltime;


        /* call setup to allocate arrays for Taylor expansions and setup global vars */
        setup_yuk(parlist, numpars, order, xyzminmax);
        thetasq = theta*theta;

        time1 = time(NULL);

        /* set global variables to track tree levels during construction */
        level = 0;
        minlevel = 50000;
        maxlevel = 0;

        printf("Creating tree... \n\n");


        create_tree_n0(&troot, 1, numpars, parlist, maxparnode,
                       xyzminmax, level, numpars);



        time2 = time(NULL);
        totaltime = difftime(time2, time1);
        *timetree = totaltime;


        printf("Tree created.\n\n");
        printf("Tree information: \n\n");

        printf("       numpar: %d\n", troot->numpar);
        printf("        x_mid: %e\n", troot->x_mid);
        printf("        y_mid: %e\n", troot->y_mid);
        printf("        z_mid: %e\n\n", troot->z_mid);
        printf("       radius: %f\n", troot->radius);
        printf("       torder: %d\n", torder);
        printf("        theta: %f\n", theta);
        printf("   maxparnode: %d\n", maxparnode);


        time1 = time(NULL);
        pc_treecode_yuk(troot, parlist, forlist,
                        numpars, kappa);
    

        time2 = time(NULL);
        totaltime = difftime(time2, time1);
    
        printf("\nTreecode time (s): %f\n", *timetree);
        printf("TOTAL time (creation + computation) (s): %f\n\n", *timetree += totaltime);

    
        /* Filling reordering array */
        for (i = 0; i < numpars; i++) orderind[orderarr[i]-1] = i;
    
    
        printf("Deallocating tree structure... \n\n");

        cleanup(troot);

        return;

} /* END function pc_treecode */


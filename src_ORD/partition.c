#include <stdio.h>

#include "tools.h"
#include "structs.h"
#include "partition.h"

/* 
 * definition of partition function
 *
 * partition determines the index MIDIND, after partitioning in place the arrays a, b, c,
 * and q, such that a(ibeg:midind) <= val and a(midind+1:iend) > val. If on entry, ibeg >
 * iend, or a(ibeg:iend) > val then midind is returned as ibeg-1.
 */

void partition(struct par *parlist, int adim, int *indarr,
               int ibeg, int iend, double val, int *midind)
{
    /* local variables */
    double ta;
    int lower, upper, tind;
    struct par temppar;


    if (ibeg < iend) {
    /*
    * temporarily stores ibeg entries and set a(ibeg) = val
    * for the partitioning algorithm.
    */
        ta = parlist[ibeg-1].r[adim];
        
        temppar = parlist[ibeg-1];
        tind = indarr[ibeg - 1];
        
        parlist[ibeg-1].r[adim] = val;
        upper = ibeg;
        lower = iend;

        while (upper != lower) {
            while ((upper < lower) && (val < parlist[lower-1].r[adim])) {
                lower--;
            }
            
            if (upper != lower) {
                parlist[upper-1] = parlist[lower-1];
                indarr[upper-1] = indarr[lower-1];
            }

            while ((upper < lower) && (val >= parlist[upper-1].r[adim])) {
                upper++;
            }

            if (upper != lower) {
                parlist[lower-1] = parlist[upper-1];
                indarr[lower-1] = indarr[upper-1];
            }
        }

        *midind = upper;


        /* replace TA in position upper and change midind if ta > val */
        if (ta > val)
            *midind = upper - 1;
        
        parlist[upper-1] = temppar;
        indarr[upper-1] = tind;

    } else if (ibeg == iend) {

        if (parlist[ibeg-1].r[adim] < val)
            *midind = ibeg;
        else
            *midind = ibeg - 1;

    } else {
        *midind = ibeg - 1;
    }

    return;

} /* END function partition */

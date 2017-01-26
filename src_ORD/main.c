#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "globtreevars.h"

#include "array.h"
#include "tools.h"
#include "routines.h"
#include "treecode.h"
#include "treecomm.h"


/* The treedriver routine in Fortran */
int main(int argc, char** argv)
{
    
    int rank, p;
    MPI_Datatype MPI_PARTYPE, MPI_TNODETYPE;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    
    
    /* runtime parameters */
    int numparsS;               /* number of particles */
    int pot_type;               /* 0 is Coulomb, 1 is screened Coulomb */
    int nt;                     /* number of time steps */
    double dt;                  /* size of time steps */
    double kappa;               /* screening parameter */
    char *sampin1, *sampout;    /* input and output files */
    
    /* treecode runtime parameters */
    int order;                  /* Taylor order for cluster interactions */
    double theta;               /* multipole acceptance criterion */
    int maxparnode;             /* maximum particles per leaf */
    
    
    /* arrays for coordinates, charges, potentials of particles */
    struct par *parlist;

    /* variables for date-time calculation */
    double t1, t2, ttemp;
    double ttree = 0.0, tcomm = 0.0, tcomp = 0.0;
    
    /* local variables */
    int i, j;
    double tpeng, tpengglob;
    int numparloc, globparloc;
    int maxparloc;
    int totexchange;
    
    int numbis;
    char binrank[20];
    double pminmax[6];
    
    struct tnode *troot;
    struct tnode *globtrees;
    MPI_Comm *commarr;
    
    
    /* Setting runtime parameters */
    sampin1 = argv[1];
    sampout = argv[2];
    numparsS = atoi(argv[3]);
    kappa = atof(argv[4]);
    pot_type = atoi(argv[5]);
    dt = atof(argv[6]);
    nt = atoi(argv[7]);
    order = atoi(argv[8]);
    theta = atof(argv[9]);
    maxparnode = atoi(argv[10]);
    
    construct_MPI_PARTYPE(&MPI_PARTYPE);
    construct_MPI_TNODETYPE(&MPI_TNODETYPE);
    
    /* Setting local variables */
    numparloc = (int)floor((double)numparsS/(double)p);
    maxparloc = numparloc + (numparsS - (int)floor((double)numparsS/(double)p) * p);
    if (rank == 0) numparloc = maxparloc;
    
    globparloc = maxparloc + numparloc * (rank-1);
    numbis = (int)log2(p);

    
    /* Allocating local structures */
    parlist = malloc(sizeof(struct par) * numparloc);
    globtrees = malloc(sizeof(struct tnode) * p);
    for (i = 0; i < p; i++) {
        make_3array((globtrees[i]).ms, order+1, order+1, order+1);
        for (j = 0; j < 8; j++) {
            globtrees[i].child[j] = malloc(sizeof(struct tnode));
            make_3array((globtrees[i].child[j])->ms, order+1, order+1, order+1);
        }
    }

    
    /* Reading in data and setting up comm network */
    readin_parlist(parlist, numparloc, globparloc, sampin1);
    contsruct_commarr(rank, numbis, binrank, &commarr);

    /* Exchanging particles to decompose particle domains */
    perform_ORB(&parlist, &numparloc, binrank, numbis, commarr, MPI_PARTYPE);
    
    
    /* Allocate and initialize global treecode data structures and variables */
    setup_yuk(parlist, numparloc, order, pminmax, theta);
    
    /* Construct tree */
    create_tree_n0(&troot, 1, numparloc, parlist, maxparnode,
                   pminmax, 0, numparloc);
    
    
    /* Receive trees from other processors */
    gather_all_trees(parlist, troot, globtrees, p, MPI_TNODETYPE);


    /* Calculate forces on local particles from local treecode and foreign trees */
    pc_treecode_yuk(troot, parlist, numparloc, kappa);
    
    pc_external_trees(globtrees, parlist, numparloc, kappa, rank, p);
    cleanup(troot);
    
    
    t1 = MPI_Wtime();
    
    /* Velocity Verlet iteration */
    for (j = 0; j < nt; j++) {
    
        for (i = 0; i < numparloc; i++) {
            parlist[i].p[0] += 0.5 * dt * parlist[i].f[0];
            parlist[i].p[1] += 0.5 * dt * parlist[i].f[1];
            parlist[i].p[2] += 0.5 * dt * parlist[i].f[2];
            
            parlist[i].r[0] += dt * parlist[i].p[0] / parlist[i].m;
            parlist[i].r[1] += dt * parlist[i].p[1] / parlist[i].m;
            parlist[i].r[2] += dt * parlist[i].p[2] / parlist[i].m;
        }
        
        
        ttemp = MPI_Wtime();
        /* Exchanging particles to decompose particle domains */
        perform_ORB(&parlist, &numparloc, binrank, numbis, commarr, MPI_PARTYPE);
        if (rank == 0) tcomm += MPI_Wtime() - ttemp;
        
        
        ttemp = MPI_Wtime();
        /* Allocate and initialize global treecode data structures and variables */
        setup_yuk(parlist, numparloc, order, pminmax, theta);
        
        /* Construct tree */
        create_tree_n0(&troot, 1, numparloc, parlist, maxparnode,
                       pminmax, 0, numparloc);
        ttree += MPI_Wtime() - ttemp;
        
        
        ttemp = MPI_Wtime();
        /* Receive trees from other processors */
        gather_all_trees(parlist, troot, globtrees, p, MPI_TNODETYPE);
        tcomm += MPI_Wtime() - ttemp;
        
        
        ttemp = MPI_Wtime();
        /* Calculate forces on local particles from local treecode and foreign trees */
        pc_treecode_yuk(troot, parlist, numparloc, kappa);
        pc_external_trees(globtrees, parlist, numparloc, kappa, rank, p);
        cleanup(troot);
        tcomp += MPI_Wtime() - ttemp;
        
        
        for (i = 0; i < numparloc; i++) {
            parlist[i].p[0] += 0.5 * dt * parlist[i].f[0];
            parlist[i].p[1] += 0.5 * dt * parlist[i].f[1];
            parlist[i].p[2] += 0.5 * dt * parlist[i].f[2];
        }
        
        tpeng = 0.0;
        for (i = 0; i < numparloc; i++) tpeng += parlist[i].peng;
        MPI_Reduce(&tpeng, &tpengglob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&exchange, &totexchange, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        
        if (rank == 0) printf("Iteration # %d... Exchanges: %d,  Potential: %e\n",
                              j, totexchange, tpengglob);
    }
    
    
    
    t2 = MPI_Wtime() - t1;
    
    //if (rank == 0) {
        printf("\n     Elapsed total time (s) on %d:  %f\n", rank, t2);
        printf(" Tree construction time (s) on %d:  %f\n", rank, ttree);
        printf("     Communication time (s) on %d:  %f\n\n", rank, tcomm);
        printf(" Force computation time (s) on %d:  %f\n\n", rank, tcomp);
    //}
    
    
    free(parlist);
    free_vector(commarr);
    
    MPI_Finalize();
    return 0;

}

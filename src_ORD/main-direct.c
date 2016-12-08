#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "array.h"
#include "tools.h"

void direct_forces(double *pS,  double *qS,  double *fS,  double *denergy,  int numpar,
                   double *pS2, double *qS2, double *fS2, double *denergy2, int numpar2,
                   int pot_type, double kappa);

void direct_forces_within(double *pS,  double *qS,  double *fS,  double *denergy,
                          int numpar, int pot_type, double kappa);

/* The treedriver routine in Fortran */
int main(int argc, char** argv)
{
    
    int rank, p;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    
    char binrank[20], binrank2[20];
    int numbis = (int)log2(p);
    int dirbis, numparglob, color;
    int numex, numex2, rank2;
    double locsum, globavg;
    double pminmax[6];
    
    struct par *parlist, *parex;
    void *tempptr;
    
    MPI_Comm commarr[numbis];
    MPI_File fp;
    MPI_Status status;
    
    /*Construct MPI datatype for particles */
    MPI_Datatype MPI_PARTYPE, oldtypes[1];
    MPI_Aint offsets[1], extent, lb;
    int blockcounts[1];
    
    offsets[0] = 0;
    oldtypes[0] = MPI_DOUBLE;
    blockcounts[0] = 8;
    MPI_Type_get_extent(MPI_DOUBLE, &lb, &extent);
    
    MPI_Type_create_struct(1, blockcounts, offsets, oldtypes, &MPI_PARTYPE);
    MPI_Type_commit(&MPI_PARTYPE);
    
    
    
    /* runtime parameters */
    int numparsS;               /* number of particles */
    int pot_type;               /* 0 is Coulomb, 1 is screened Coulomb */
    int nt;                     /* number of time steps */
    double dt;                  /* size of time steps */
    double kappa;               /* screening parameter */
    char *sampin1, *sampout;    /* input and output files */
    
    /* arrays for coordinates, charges, potentials of particles */
    struct foreng *forenglist;

    
    /* variables for date-time calculation */
    double t1, t2;
    
    /* local variables */
    int i, j, k;
    int xdim;
    double dengyglob, dengyloc;
    double buf[5];
    int numparloc, globparloc;
    int maxparloc, senditer;
    int numpar1, numpar2, sendtag;
    
    
    sampin1 = argv[1];
    sampout = argv[2];
    numparsS = atoi(argv[3]);
    kappa = atof(argv[4]);
    pot_type = atoi(argv[5]);
    dt = atof(argv[6]);
    nt = atoi(argv[7]);

    numparloc = (int)floor((double)numparsS/(double)p);
    maxparloc = numparloc + (numparsS - (int)floor((double)numparsS/(double)p) * p);
    xdim = maxparloc * 3;
    
    if (rank == 0)
        numparloc = maxparloc;
    
    globparloc = maxparloc + numparloc * (rank-1);
    
    
    printf("Proc %d has %d local particles\n", rank, numparloc);
    
    
    /* Allocate arrays */
    /* Allocate all maxparloc for sending equality */
    make_vector(parlist, maxparloc);
    make_vector(parex, maxparloc);
    

    MPI_File_open(MPI_COMM_WORLD, sampin1, MPI_MODE_RDONLY, MPI_INFO_NULL, &fp);
    MPI_File_seek(fp, (MPI_Offset)globparloc*5*sizeof(double), MPI_SEEK_SET);
    
    for (i = 0; i < numparloc; i++) {
        MPI_File_read(fp, buf, 5, MPI_DOUBLE, &status);
        parlist[i].r[0] = buf[0];
        parlist[i].r[1] = buf[1];
        parlist[i].r[2] = buf[2];
        parlist[i].q = buf[3];
        parlist[i].m = buf[4];
        parlist[i].p[0] = 0.0;
        parlist[i].p[1] = 0.0;
        parlist[i].p[2] = 0.0;
    }
    
    MPI_File_close(&fp);

    
    
    

    
    inttobinstr(rank, numbis, binrank);

    commarr[0] = MPI_COMM_WORLD;
    for (i=0; i < numbis-1; i++) {
        if (binrank[i] == '0') color = 0;
        else color = 1;
        
        MPI_Comm_split(commarr[i], color, rank, &commarr[i+1]);
    }
    
    /*
    for (i=0; i < numbis; i++) {
        MPI_Comm_size(commarr[i], &color);
        printf("On proc %d, comm %d has size: %d\n", rank, i, color);
    }
    */
    
    
    for (i=0; i < numbis; i++) {
        dirbis = i%3;
        
        locsum = 0;
        for (j=0; j < numparloc; j++) {
            locsum += parlist[j].r[dirbis];
        }
        
        MPI_Allreduce(&locsum, &globavg, 1, MPI_DOUBLE, MPI_SUM, commarr[i]);
        MPI_Allreduce(&numparloc, &numparglob, 1, MPI_INT, MPI_SUM, commarr[i]);
        globavg /= (double)numparglob;
        

        /*
        if (binrank[i] == '0')
            printf("Proc %d, binrank %s must be below %f in dim %d\n", rank, binrank, globavg, i);
        else if (binrank[i] == '1')
            printf("Proc %d, binrank %s must be above %f in dim %d\n", rank, binrank, globavg, i);
        }
        */
   
   
        
        if (binrank[i] == '0') {
            numex = 0;
            for (j = 0; j < numparloc; j++) {
                if (parlist[j].r[dirbis] > globavg) {
                    parex[numex] = parlist[j];
                    remove_particle(parlist, j, numparloc);
                    
                    numex++; j--; numparloc--;
                }
            }
        } else {
            numex = 0;
            for (j = 0; j < numparloc; j++) {
                if (parlist[j].r[dirbis] <= globavg) {
                    parex[numex] = parlist[j];
                    remove_particle(parlist, j, numparloc);
                    
                    numex++; j--; numparloc--;
                }
            }
        }
        

        /*
        for (j=0; j < numparloc; j++) {
            printf("    Left in: %d, %s: %d %f %f %f %f\n", rank, binrank, j, parlist[j].r[0], parlist[j].r[1],
                   parlist[j].r[2], parlist[j].m);
        }
        
        for (j=0; j < numex; j++) {
            printf("    Ship on: %d, %s: %d %f %f %f %f\n", rank, binrank, j, parex[j].r[0], parex[j].r[1],
                   parex[j].r[2], parex[j].m);
        }
        */
        
      
        
        strncpy(binrank2, binrank, 20);
        if (binrank2[i] == '0') binrank2[i] = '1';
        else if (binrank2[i] == '1') binrank2[i] = '0';
        
        rank2 = strtol(binrank2,NULL,2);
        /*
        printf("Proc %d, %s: %d particles staged for exchange to %d, %s across %d, %d left local\n",
               rank, binrank, numex, rank2, binrank2, i, numparloc);
        */
        
        
        MPI_Sendrecv(&numex,  1, MPI_INT, rank2, 0,
                     &numex2, 1, MPI_INT, rank2, 0,
                     MPI_COMM_WORLD, &status);
        /*
        printf("Proc %d, %s: proc %d, %s wants to send %d particles over\n",
               rank, binrank, rank2, binrank2, numex2);
        */
        
        tempptr = realloc(parlist, (numparloc+numex2)*sizeof(struct par));
        parlist = tempptr;
        
        MPI_Sendrecv(parex, numex, MPI_PARTYPE, rank2, 0,
                     &parlist[numparloc], numex2, MPI_PARTYPE, rank2, 0,
                     MPI_COMM_WORLD, &status);
        
        numparloc += numex2;
        
        /*
        for (j = 0; j < numparloc; j++)
            printf("    Proc %d, %s, after bis %d: %d %f %f %f %f %p\n", rank, binrank, i, j,
                   parlist[j].r[0], parlist[j].r[1],
                   parlist[j].r[2], parlist[j].m, &parlist[j]);
        */

    }

    
    
    for (i = 0; i < 3; i++) pminmax[2*i]   = parlist[0].r[i];
    for (i = 0; i < 3; i++) pminmax[2*i+1] = parlist[0].r[i];
    
    for (i = 0; i < numparloc; i++) {
        pminmax[0] = fmin(pminmax[0], parlist[i].r[0]);
        pminmax[1] = fmax(pminmax[1], parlist[i].r[0]);
        pminmax[2] = fmin(pminmax[2], parlist[i].r[1]);
        pminmax[3] = fmax(pminmax[3], parlist[i].r[1]);
        pminmax[4] = fmin(pminmax[4], parlist[i].r[2]);
        pminmax[5] = fmax(pminmax[5], parlist[i].r[2]);
    }

    printf("Proc %d now has %d particles, inside box: %f, %f, %f, %f, %f, %f.\n",
           rank, numparloc, minx, maxx, miny, maxy, minz, maxz);
    
    
    
    
    
    
    MPI_Finalize();
    return 0;

}




void direct_forces(double *pS,  double *qS,  double *fS,  double *denergy,  int numpar,
                   double *pS2, double *qS2, double *fS2, double *denergy2, int numpar2,
                   int pot_type, double kappa)
{
    /* local variables */
    int i, j;
    double tx, ty, tz, xi, yi, zi, rad;
    double temp, temp1, temp2;
    
    
    if (pot_type == 0) {
        
        for (j = 0; j < numpar2; j++) {
            fS2[j*3 + 0] = 0.0;  fS2[j*3 + 1] = 0.0;  fS2[j*3 + 2] = 0.0;
            denergy2[j] = 0.0;
        }

        for (i = 0; i < numpar; i++) {
            xi = pS[i*3 + 0];
            yi = pS[i*3 + 1];
            zi = pS[i*3 + 2];
            
            for (j = 0; j < numpar2; j++) {
                tx = xi - pS2[j*3 + 0];
                ty = yi - pS2[j*3 + 1];
                tz = zi - pS2[j*3 + 2];

                rad = sqrt(tx*tx + ty*ty + tz*tz);
                    
                denergy[i] = denergy[i] + qS2[j] / rad;
                denergy2[j] = denergy2[j] + qS[i] / rad;

                temp1 = qS2[j] / pow(rad, 3);
                temp2 = qS[i] / pow(rad, 3);
                
                fS[i*3 + 0] += tx * temp1;
                fS[i*3 + 1] += ty * temp1;
                fS[i*3 + 2] += tz * temp1;

                fS2[j*3 + 0] -= tx * temp2;
                fS2[j*3 + 1] -= ty * temp2;
                fS2[j*3 + 2] -= tz * temp2;
            }
        }
        
    } else if (pot_type == 1) {
        
        for (j = 0; j < numpar2; j++) {
            fS2[j*3 + 0] = 0.0;  fS2[j*3 + 1] = 0.0;  fS2[j*3 + 2] = 0.0;
            denergy2[j] = 0.0;
        }

        for (i = 0; i < numpar; i++) {
            xi = pS[i*3 + 0];
            yi = pS[i*3 + 1];
            zi = pS[i*3 + 2];
            
            
            for (j = 0; j < numpar2; j++) {
                tx = xi - pS2[j*3 + 0];
                ty = yi - pS2[j*3 + 1];
                tz = zi - pS2[j*3 + 2];

                rad = sqrt(tx*tx + ty*ty + tz*tz);
                temp = exp(-kappa * rad);
                
                denergy[i] = denergy[i] + qS2[j] * temp / rad;
                denergy2[j] = denergy2[j] + qS[i] * temp / rad;

                temp1 = qS2[j] * temp * (kappa * rad + 1) / pow(rad, 3);
                temp2 = qS[i] * temp * (kappa * rad + 1) / pow(rad, 3);

                fS[i*3 + 0] += temp1 * tx;
                fS[i*3 + 1] += temp1 * ty;
                fS[i*3 + 2] += temp1 * tz;

                fS2[j*3 + 0] -= temp2 * tx;
                fS2[j*3 + 1] -= temp2 * ty;
                fS2[j*3 + 2] -= temp2 * tz;
            }
        }
    }

    return;
}



void direct_forces_within(double *pS,  double *qS,  double *fS,  double *denergy,
                          int numpar, int pot_type, double kappa)
{
    /* local variables */
    int i, j;
    double tx, ty, tz, xi, yi, zi, rad;
    double temp, temp1, temp2;
    
    if (pot_type == 0) {
        
        for (i = 0; i < numpar; i++) {
            fS[i*3 + 0] = 0.0;  fS[i*3 + 1] = 0.0;  fS[i*3 + 2] = 0.0;
            denergy[i] = 0.0;
        }
        
        for (i = 0; i < numpar; i++) {
            xi = pS[i*3 + 0];
            yi = pS[i*3 + 1];
            zi = pS[i*3 + 2];
            
            for (j = i+1; j < numpar; j++) {
                tx = xi - pS[j*3 + 0];
                ty = yi - pS[j*3 + 1];
                tz = zi - pS[j*3 + 2];
                
                rad = sqrt(tx*tx + ty*ty + tz*tz);
                
                denergy[i] = denergy[i] + qS[j] / rad;
                denergy[j] = denergy[j] + qS[i] / rad;
                
                temp1 = qS[j] / pow(rad, 3);
                temp2 = qS[i] / pow(rad, 3);
                
                fS[i*3 + 0] += tx * temp1;
                fS[i*3 + 1] += ty * temp1;
                fS[i*3 + 2] += tz * temp1;
                
                fS[j*3 + 0] -= tx * temp2;
                fS[j*3 + 1] -= ty * temp2;
                fS[j*3 + 2] -= tz * temp2;
            }
        }
        
    } else if (pot_type == 1) {
        
        for (i = 0; i < numpar; i++) {
            fS[i*3 + 0] = 0.0;  fS[i*3 + 1] = 0.0;  fS[i*3 + 2] = 0.0;
            denergy[i] = 0.0;
        }
        
        for (i = 0; i < numpar; i++) {
            xi = pS[i*3 + 0];
            yi = pS[i*3 + 1];
            zi = pS[i*3 + 2];
            
            for (j = i+1; j < numpar; j++) {
                tx = xi - pS[j*3 + 0];
                ty = yi - pS[j*3 + 1];
                tz = zi - pS[j*3 + 2];
                
                rad = sqrt(tx*tx + ty*ty + tz*tz);
                temp = exp(-kappa * rad);
                
                denergy[i] = denergy[i] + qS[j] * temp / rad;
                denergy[j] = denergy[j] + qS[i] * temp / rad;
                
                temp1 = qS[j] * temp * (kappa * rad + 1) / pow(rad, 3);
                temp2 = qS[i] * temp * (kappa * rad + 1) / pow(rad, 3);
                
                fS[i*3 + 0] += temp1 * tx;
                fS[i*3 + 1] += temp1 * ty;
                fS[i*3 + 2] += temp1 * tz;
                
                fS[j*3 + 0] -= temp2 * tx;
                fS[j*3 + 1] -= temp2 * ty;
                fS[j*3 + 2] -= temp2 * tz;
            }
        }
    }
    
    return;
}

#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
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
    
    /* runtime parameters */
    int numparsS;               /* number of particles */
    int pot_type;               /* 0 is Coulomb, 1 is screened Coulomb */
    int nt;                     /* number of time steps */
    double dt;                  /* size of time steps */
    double kappa;               /* screening parameter */
    char *sampin1, *sampout;    /* input and output files */
    
    /* arrays for coordinates, charges, potentials of particles */
    double *qS, *mS;            /* particle charges and masses */
    double *rS;               /* particle positions */
    double *fS;               /* particle forces */
    double *pS;                /* particle momenta */
    double *denergy;           /* particle pot energies */
    double *kenergy;           /* particle kin energies */
    
    double *qS2, *mS2;          /* particle charges and masses */
    double *rS2;              /* particle positions */
    double *fS2;              /* particle forces */
    double *denergy2;          /* particle pot energies */
    
    double *rS1, *rStemp;
    double *fS1;
    double *qS1, *qStemp;
    double *denergy1;
    
    /* variables for date-time calculation */
    double t1, t2, ttemp;
    double tcomm = 0.0;
    
    /* local variables */
    int i, j, k;
    int xdim;
    double dengyglob, dengyloc;
    double buf[5];
    int numparloc;
    int globparloc;
    int maxparloc;
    int senditer;

    int numpar1;
    int numpar2;
    int sendtag;
    
    /* MPI Variables */
    MPI_File fp;
    MPI_Status status;
    
    
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
    
    if (p%2 == 1)
        senditer = (int)((p-1)/2);          /* for odd p, all procs do senditer */
    else
        senditer = (int)(floor((p-1)/2));   /*for even p, half do senditer, half 1 more */
    
    
    printf("Proc %d has %d local particles\n", rank, numparloc);
    
    
    /* Allocate arrays */
    /* Allocate all maxparloc for sending equality */
    
    make_vector(rS, (nt+1) * maxparloc * 3);
    make_vector(fS, (nt+1) * maxparloc * 3);
    
    make_vector(denergy, (nt+1) * maxparloc);
    make_vector(kenergy, (nt+1) * maxparloc);
    
    make_vector(rS2, maxparloc * 3);
    make_vector(fS2, maxparloc * 3);
    make_vector(denergy2, maxparloc);
    
    make_vector(rS1, maxparloc * 3);
    make_vector(fS1, maxparloc * 3);
    make_vector(denergy1, maxparloc);
    
    make_vector(qS, maxparloc);
    make_vector(mS, maxparloc);
    make_vector(pS, maxparloc * 3);
    
    make_vector(qS1, maxparloc);
    make_vector(qS2, maxparloc);
    make_vector(mS2, maxparloc);
    
    
    MPI_File_open(MPI_COMM_WORLD, sampin1, MPI_MODE_RDONLY, MPI_INFO_NULL, &fp);
    MPI_File_seek(fp, (MPI_Offset)globparloc*5*sizeof(double), MPI_SEEK_SET);
    
    for (i = 0; i < numparloc; i++) {
        MPI_File_read(fp, buf, 5, MPI_DOUBLE, &status);
        rS[i*3 + 0] = buf[0];
        rS[i*3 + 1] = buf[1];
        rS[i*3 + 2] = buf[2];
        qS[i] = buf[3];
        mS[i] = buf[4];
    }
    
    MPI_File_close(&fp);


    /* Zero out momenta before dynamics simulation */
    for (i = 0; i < numparloc; i++) {
        pS[i*3 + 0] = 0.0;  pS[i*3 + 1] = 0.0;  pS[i*3 + 2] = 0.0;
    }
    

    /* Calculation of forces at initial time */
    direct_forces_within(&rS[0], qS, &fS[0], &denergy[0],
                         numparloc, pot_type, kappa);
    
    for (i = 0; i < numparloc; i++) {
        rS1[i*3 + 0] = rS[0 + i*3 + 0];
        rS1[i*3 + 1] = rS[0 + i*3 + 1];
        rS1[i*3 + 2] = rS[0 + i*3 + 2];
        qS1[i] = qS[i];
    }

    numpar1 = numparloc;
    sendtag = rank;
    
    for (k = 0; k < senditer; k++) {
        MPI_Sendrecv(rS1, xdim, MPI_DOUBLE, mod(rank+1,p), sendtag,      /*stag is origin*/
                     rS2, xdim, MPI_DOUBLE, mod(rank-1,p), MPI_ANY_TAG,  /*rtag is origin*/
                     MPI_COMM_WORLD, &status);
        
        MPI_Sendrecv(qS1, maxparloc, MPI_DOUBLE, mod(rank+1,p), sendtag,
                     qS2, maxparloc, MPI_DOUBLE, mod(rank-1,p), MPI_ANY_TAG,
                     MPI_COMM_WORLD, &status);
        
        MPI_Sendrecv(&numpar1, 1, MPI_INT, mod(rank+1,p), sendtag,
                     &numpar2, 1, MPI_INT, mod(rank-1,p), MPI_ANY_TAG,
                     MPI_COMM_WORLD, &status);
        
        sendtag = status.MPI_TAG;
        
        direct_forces(&rS[0],  qS,  &fS[0], &denergy[0], numparloc,
                         rS2, qS2,     fS2,    denergy2, numpar2,
                      pot_type, kappa);
        
        MPI_Sendrecv(fS2, xdim, MPI_DOUBLE, sendtag, sendtag,        /*stag is origin*/
                     fS1, xdim, MPI_DOUBLE, mod(rank+k+1,p), rank,   /*rtag is origin*/
                     MPI_COMM_WORLD, &status);
        
        MPI_Sendrecv(denergy2, maxparloc, MPI_DOUBLE, sendtag, sendtag,        /*stag is origin*/
                     denergy1, maxparloc, MPI_DOUBLE, mod(rank+k+1,p), rank,   /*rtag is origin*/
                     MPI_COMM_WORLD, &status);
        
        for (i = 0; i < numparloc; i++) {
            fS[0 + i*3 + 0] += fS1[i*3 + 0];
            fS[0 + i*3 + 1] += fS1[i*3 + 1];
            fS[0 + i*3 + 2] += fS1[i*3 + 2];
            denergy[0 + i] += denergy1[i];
        }
        
        rStemp = rS1; qStemp = qS1;
        rS1 = rS2;    qS1 = qS2;
        rS2 = rStemp; qS2 = qStemp;
        
        numpar1 = numpar2;
    }
    
    /* For even number of procs, must perform a half iteration */
    if (p%2 == 0) {
        if(rank >= senditer && rank < p/2+senditer) {
            MPI_Send(rS1, xdim, MPI_DOUBLE, mod(rank+1,p), sendtag, MPI_COMM_WORLD);
            MPI_Send(qS1, maxparloc, MPI_DOUBLE, mod(rank+1,p), sendtag, MPI_COMM_WORLD);
            MPI_Send(&numpar1, 1, MPI_INT, mod(rank+1,p), sendtag, MPI_COMM_WORLD);
        }
        
        if(rank >= p/2) {
            MPI_Recv(rS2, xdim, MPI_DOUBLE, mod(rank-1,p), MPI_ANY_TAG,
                     MPI_COMM_WORLD, &status);
            MPI_Recv(qS2, maxparloc, MPI_DOUBLE, mod(rank-1,p), MPI_ANY_TAG,
                     MPI_COMM_WORLD, &status);
            MPI_Recv(&numpar2, 1, MPI_INT, mod(rank-1,p), MPI_ANY_TAG,
                     MPI_COMM_WORLD, &status);
            
            sendtag = status.MPI_TAG;
            
            direct_forces(&rS[0],  qS,  &fS[0], &denergy[0], numparloc,
                          rS2, qS2,     fS2,    denergy2, numpar2,
                          pot_type, kappa);
            
            MPI_Send(fS2, xdim, MPI_DOUBLE, sendtag, sendtag, MPI_COMM_WORLD);
            MPI_Send(denergy2, maxparloc, MPI_DOUBLE, sendtag, sendtag, MPI_COMM_WORLD);
        }
        
        if(rank < p/2) {
            MPI_Recv(fS1, xdim, MPI_DOUBLE, mod(rank+senditer+1,p), rank,
                     MPI_COMM_WORLD, &status);
            MPI_Recv(denergy1, maxparloc, MPI_DOUBLE, mod(rank+senditer+1,p), rank,
                     MPI_COMM_WORLD, &status);
            
            for (i = 0; i < numparloc; i++) {
                fS[0 + i*3 + 0] += fS1[i*3 + 0];
                fS[0 + i*3 + 1] += fS1[i*3 + 1];
                fS[0 + i*3 + 2] += fS1[i*3 + 2];
                denergy[0 + i] += denergy1[i];
            }
        }
    }
    
    for (i = 0; i < numparloc; i++) {
        fS[0 + i*3 + 0] *= qS[i];
        fS[0 + i*3 + 1] *= qS[i];
        fS[0 + i*3 + 2] *= qS[i];
        denergy[0 + i] *= qS[i];
    }


    /* Starting timer before dynamics simulation */
    if (rank == 0)
        t1 = MPI_Wtime();

    
    /* Velocity Verlet iteration */
    for (j = 0; j < nt; j++) {
        
        for (i = 0; i < numparloc; i++) {

            pS[i*3+0] += 0.5 * dt * fS[j*xdim + i*3 + 0];
            pS[i*3+1] += 0.5 * dt * fS[j*xdim + i*3 + 1];
            pS[i*3+2] += 0.5 * dt * fS[j*xdim + i*3 + 2];
            
            rS[(j+1)*xdim + i*3 + 0] = rS[j*xdim + i*3 + 0] + dt * pS[i*3+0] / mS[i];
            rS[(j+1)*xdim + i*3 + 1] = rS[j*xdim + i*3 + 1] + dt * pS[i*3+1] / mS[i];
            rS[(j+1)*xdim + i*3 + 2] = rS[j*xdim + i*3 + 2] + dt * pS[i*3+2] / mS[i];
        }
        
        
 
        
        direct_forces_within(&rS[(j+1)*(maxparloc*3)], qS,
                             &fS[(j+1)*(maxparloc*3)],
                             &denergy[(j+1)*maxparloc],
                             numparloc, pot_type, kappa);
        
        for (i = 0; i < numparloc; i++) {
            rS1[i*3+0] = rS[(j+1)*xdim + i*3 + 0];
            rS1[i*3+1] = rS[(j+1)*xdim + i*3 + 1];
            rS1[i*3+2] = rS[(j+1)*xdim + i*3 + 2];
            qS1[i] = qS[i];
        }
        
        numpar1 = numparloc;
        sendtag = rank;
        
        for (k = 0; k < senditer; k++) {
            
            if (rank == 0) ttemp = MPI_Wtime();
            MPI_Sendrecv(rS1, xdim, MPI_DOUBLE, mod(rank+1,p), sendtag,      /*stag is origin*/
                         rS2, xdim, MPI_DOUBLE, mod(rank-1,p), MPI_ANY_TAG,  /*rtag is origin*/
                         MPI_COMM_WORLD, &status);
            
            MPI_Sendrecv(qS1, maxparloc, MPI_DOUBLE, mod(rank+1,p), sendtag,
                         qS2, maxparloc, MPI_DOUBLE, mod(rank-1,p), MPI_ANY_TAG,
                         MPI_COMM_WORLD, &status);
            
            MPI_Sendrecv(&numpar1, 1, MPI_INT, mod(rank+1,p), sendtag,
                         &numpar2, 1, MPI_INT, mod(rank-1,p), MPI_ANY_TAG,
                         MPI_COMM_WORLD, &status);
            
            sendtag = status.MPI_TAG;
            if (rank == 0) tcomm += MPI_Wtime() - ttemp;
            
            
            direct_forces(&rS[(j+1)*xdim], qS, &fS[(j+1)*xdim],
                          &denergy[(j+1)*maxparloc], numparloc,
                          rS2, qS2, fS2, denergy2, numpar2,
                          pot_type, kappa);

            
            if (rank == 0) ttemp = MPI_Wtime();
            MPI_Sendrecv(fS2, xdim, MPI_DOUBLE, sendtag, sendtag,        /*stag is origin*/
                         fS1, xdim, MPI_DOUBLE, mod(rank+k+1,p), rank,   /*rtag is origin*/
                         MPI_COMM_WORLD, &status);
            
            MPI_Sendrecv(denergy2, maxparloc, MPI_DOUBLE, sendtag, sendtag,        /*stag is origin*/
                         denergy1, maxparloc, MPI_DOUBLE, mod(rank+k+1,p), rank,   /*rtag is origin*/
                         MPI_COMM_WORLD, &status);
            if (rank == 0) tcomm += MPI_Wtime() - ttemp;
            
            
            for (i = 0; i < numparloc; i++) {
                fS[(j+1)*xdim + i*3 + 0] += fS1[i*3 + 0];
                fS[(j+1)*xdim + i*3 + 1] += fS1[i*3 + 1];
                fS[(j+1)*xdim + i*3 + 2] += fS1[i*3 + 2];
                denergy[(j+1)*maxparloc + i] += denergy1[i];
            }
            
            rStemp = rS1; qStemp = qS1;
            rS1 = rS2;    qS1 = qS2;
            rS2 = rStemp; qS2 = qStemp;
            
            numpar1 = numpar2;
        }
        
        /* For even number of procs, must perform a half iteration */
        if (p%2 == 0) {
            
            if (rank == 0) ttemp = MPI_Wtime();
            if(rank >= senditer && rank < p/2+senditer) {
                MPI_Send(rS1, xdim, MPI_DOUBLE, mod(rank+1,p), sendtag, MPI_COMM_WORLD);
                MPI_Send(qS1, maxparloc, MPI_DOUBLE, mod(rank+1,p), sendtag, MPI_COMM_WORLD);
                MPI_Send(&numpar1, 1, MPI_INT, mod(rank+1,p), sendtag, MPI_COMM_WORLD);
            }
            if (rank == 0) tcomm += MPI_Wtime() - ttemp;
            
            
            if(rank >= p/2) {
                MPI_Recv(rS2, xdim, MPI_DOUBLE, mod(rank-1,p), MPI_ANY_TAG,
                         MPI_COMM_WORLD, &status);
                MPI_Recv(qS2, maxparloc, MPI_DOUBLE, mod(rank-1,p), MPI_ANY_TAG,
                         MPI_COMM_WORLD, &status);
                MPI_Recv(&numpar2, 1, MPI_INT, mod(rank-1,p), MPI_ANY_TAG,
                         MPI_COMM_WORLD, &status);
                
                sendtag = status.MPI_TAG;
                
                direct_forces(&rS[(j+1)*xdim], qS, &fS[(j+1)*xdim],
                              &denergy[(j+1)*maxparloc], numparloc,
                              rS2, qS2, fS2, denergy2, numpar2,
                              pot_type, kappa);
                
                MPI_Send(fS2, xdim, MPI_DOUBLE, sendtag, sendtag, MPI_COMM_WORLD);
                MPI_Send(denergy2, maxparloc, MPI_DOUBLE, sendtag, sendtag, MPI_COMM_WORLD);
            }
            
            if(rank < p/2) {
                
                if (rank == 0) ttemp = MPI_Wtime();
                MPI_Recv(fS1, xdim, MPI_DOUBLE, mod(rank+senditer+1,p), rank,
                         MPI_COMM_WORLD, &status);
                MPI_Recv(denergy1, maxparloc, MPI_DOUBLE, mod(rank+senditer+1,p), rank,
                         MPI_COMM_WORLD, &status);
                if (rank == 0) tcomm += MPI_Wtime() - ttemp;
                
                
                for (i = 0; i < numparloc; i++) {
                    fS[(j+1)*xdim + i*3 + 0] += fS1[i*3 + 0];
                    fS[(j+1)*xdim + i*3 + 1] += fS1[i*3 + 1];
                    fS[(j+1)*xdim + i*3 + 2] += fS1[i*3 + 2];
                    denergy[(j+1)*maxparloc + i] += denergy1[i];
                }
            }
        }
        
        for (i = 0; i < numparloc; i++) {
            fS[(j+1)*xdim + i*3 + 0] *= qS[i];
            fS[(j+1)*xdim + i*3 + 1] *= qS[i];
            fS[(j+1)*xdim + i*3 + 2] *= qS[i];
            denergy[(j+1)*maxparloc + i] *= qS[i];
        }

        
        
        for (i = 0; i < numparloc; i++) {
            pS[i*3+0] += 0.5 * dt * fS[(j+1)*xdim + i*3 + 0];
            pS[i*3+1] += 0.5 * dt * fS[(j+1)*xdim + i*3 + 1];
            pS[i*3+2] += 0.5 * dt * fS[(j+1)*xdim + i*3 + 2];
        }
        
        dengyloc = sum(&denergy[(j+1)*maxparloc], numparloc);
        MPI_Reduce(&dengyloc, &dengyglob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        
        if (rank == 0) printf("Iteration # %d... Global potential: %e\n", j, dengyglob);
    }

    if (rank == 0) t2 = MPI_Wtime() - t1;
    
    dengyloc = sum(&denergy[nt*maxparloc], numparloc);
    MPI_Reduce(&dengyloc, &dengyglob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        printf("\n  Elapsed total time (s): %f s\n", t2);
        printf("  Communication time (s): %f s\n", tcomm);
        printf("        Global potential: %e\n\n", dengyglob);
    }

    
    /* Deallocate arrays */
    free_vector(rS);
    free_vector(fS);
    free_vector(denergy);
    free_vector(kenergy);
    
    free_vector(qS);
    free_vector(mS);
    free_vector(pS);
    
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

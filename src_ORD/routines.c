#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "array.h"
#include "tools.h"
#include "routines.h"

/* Declare global exchange variable */
int exchange;


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




void direct_forces_within(struct par *parlist, struct foreng *forlist,
                          int numpars, int pot_type, double kappa)
{
    /* local variables */
    int i, j;
    double tx, ty, tz, xi, yi, zi, rad;
    double temp, temp1, temp2, qi;
    
    if (pot_type == 0) {
        for (i = 0; i < numpars; i++) {
            forlist[i].peng = 0.0;
            forlist[i].f[0] = 0.0;
            forlist[i].f[1] = 0.0;
            forlist[i].f[2] = 0.0;
        }
        
        for (i = 0; i < numpars; i++) {
            xi = parlist[i].r[0];
            yi = parlist[i].r[1];
            zi = parlist[i].r[2];
            qi = parlist[i].q;
            
            for (j = i+1; j < numpars; j++) {
                tx = xi - parlist[j].r[0];
                ty = yi - parlist[j].r[1];
                tz = zi - parlist[j].r[2];
                
                rad = sqrt(tx*tx + ty*ty + tz*tz);
                
                forlist[i].peng += parlist[j].q / rad;
                forlist[j].peng +=           qi / rad;
                
                temp1 = parlist[j].q / pow(rad, 3);
                temp2 =           qi / pow(rad, 3);
                
                forlist[i].f[0] += tx * temp1;
                forlist[i].f[1] += ty * temp1;
                forlist[i].f[2] += tz * temp1;
                
                forlist[j].f[0] -= tx * temp2;
                forlist[j].f[1] -= ty * temp2;
                forlist[j].f[2] -= tz * temp2;
            }
        }
        
    } else if (pot_type == 1) {
        for (i = 0; i < numpars; i++) {
            forlist[i].peng = 0.0;
            forlist[i].f[0] = 0.0;
            forlist[i].f[1] = 0.0;
            forlist[i].f[2] = 0.0;
        }
        
        for (i = 0; i < numpars; i++) {
            xi = parlist[i].r[0];
            yi = parlist[i].r[1];
            zi = parlist[i].r[2];
            qi = parlist[i].q;
            
            for (j = i+1; j < numpars; j++) {
                tx = xi - parlist[j].r[0];
                ty = yi - parlist[j].r[1];
                tz = zi - parlist[j].r[2];
                
                rad = sqrt(tx*tx + ty*ty + tz*tz);
                temp = exp(-kappa * rad);
                
                forlist[i].peng += parlist[j].q * temp / rad;
                forlist[j].peng +=           qi * temp / rad;
                
                temp1 = parlist[j].q * temp * (kappa * rad + 1) / pow(rad, 3);
                temp2 =           qi * temp * (kappa * rad + 1) / pow(rad, 3);
                
                forlist[i].f[0] += tx * temp1;
                forlist[i].f[1] += ty * temp1;
                forlist[i].f[2] += tz * temp1;
                
                forlist[j].f[0] -= tx * temp2;
                forlist[j].f[1] -= ty * temp2;
                forlist[j].f[2] -= tz * temp2;
            }
        }
    }
    
    return;
}




void perform_ORB(struct par **parlist, int *numparloc, char binrank[20], int numbis,
                 MPI_Comm commarr[numbis], MPI_Datatype MPI_PARTYPE)
{
    int i, j;
    
    char binrank2[20];
    int dirbis, numparglob;
    int numex, numex2, rank2;
    double locsum, globavg;
    
    struct par *parex;
    struct par *tempptr;
    MPI_Status status;
    
    parex = malloc(sizeof(struct par)*(*numparloc));
    
    /* Zero out global exchange tracking variable */
    exchange = 0;
    
    
    for (i=0; i < numbis; i++) {
        dirbis = i%3;
        
        locsum = 0;
        for (j=0; j < *numparloc; j++) {
            locsum += (*parlist)[j].r[dirbis];
        }
        
        MPI_Allreduce(&locsum, &globavg, 1, MPI_DOUBLE, MPI_SUM, commarr[i]);
        MPI_Allreduce(numparloc, &numparglob, 1, MPI_INT, MPI_SUM, commarr[i]);
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
            for (j = 0; j < *numparloc; j++) {
                if ((*parlist)[j].r[dirbis] > globavg) {
                    parex[numex] = (*parlist)[j];
                    remove_particle(*parlist, j, *numparloc);
                    
                    numex++; j--; (*numparloc)--;
                }
            }
        } else {
            numex = 0;
            for (j = 0; j < *numparloc; j++) {
                if ((*parlist)[j].r[dirbis] <= globavg) {
                    parex[numex] = (*parlist)[j];
                    remove_particle(*parlist, j, *numparloc);
                    
                    numex++; j--; (*numparloc)--;
                }
            }
        }
        
        /* Global variable to track # of exchanges */
        exchange += numex;
        
        
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
        
         //printf("Proc %s: proc %d, %s wants to send %d particles over\n",
         //binrank, rank2, binrank2, numex2);
        
        
        tempptr = realloc(*parlist, ((*numparloc)+numex2)*sizeof(struct par));
        *parlist = tempptr;
        
        /*
        printf("Proc %s: reallocated parlist to size %d from size %d\n", binrank,
               (*numparloc)+numex2, *numparloc);
         */
        
        
        MPI_Sendrecv(parex, numex, MPI_PARTYPE, rank2, 0,
                     &((*parlist)[*numparloc]), numex2, MPI_PARTYPE, rank2, 0,
                     MPI_COMM_WORLD, &status);
        
        (*numparloc) += numex2;
        
        
        /*
         for (j = *numparloc-10; j < *numparloc; j++)
         printf("    Proc %s, after bis %d: %d %f %f %f %f %p\n", binrank, i, j,
         (*parlist)[j].r[0], (*parlist)[j].r[1],
         (*parlist)[j].r[2], (*parlist)[j].m, &((*parlist)[j]));
         */
         
        
        
    }
    
    free(parex);
    
    return;
}




void contsruct_commarr(int rank, int numbis, char binrank[20], MPI_Comm **commarr)
{
    int i, color;
    inttobinstr(rank, numbis, binrank);
    
    make_vector(*commarr, numbis);
    (*commarr)[0] = MPI_COMM_WORLD;
    
    for (i=0; i < numbis-1; i++) {
        if (binrank[i] == '0') color = 0;
        else color = 1;
        
        MPI_Comm_split((*commarr)[i], color, rank, &(*commarr)[i+1]);
    }
    
    return;
}




void readin_parlist(struct par *parlist, int numparloc, int globparloc, char *sampin1)
{
    int i;
    double buf[5];
    MPI_File fp;
    MPI_Status status;
    
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
        
        parlist[i].f[0] = 0.0;
        parlist[i].f[1] = 0.0;
        parlist[i].f[2] = 0.0;
        parlist[i].peng = 0.0;
    }
    
    MPI_File_close(&fp);
    
    return;
}




void construct_MPI_PARTYPE(MPI_Datatype *MPI_PARTYPE)
{
    MPI_Datatype oldtypes[1];
    MPI_Aint offsets[1], extent, lb;
    int blockcounts[1];
    
    offsets[0] = 0;
    oldtypes[0] = MPI_DOUBLE;
    blockcounts[0] = 12;
    MPI_Type_get_extent(MPI_DOUBLE, &lb, &extent);
    
    MPI_Type_create_struct(1, blockcounts, offsets, oldtypes, MPI_PARTYPE);
    MPI_Type_commit(MPI_PARTYPE);
    
    return;
}




void construct_MPI_TNODETYPE(MPI_Datatype *MPI_TNODETYPE)
{
    MPI_Datatype oldtypes[3] = {MPI_INT, MPI_DOUBLE, MPI_INT};
    int blockcounts[3] = {3, 12, 4};
    MPI_Aint offsets[3], extent[2], lb;
    
    MPI_Type_get_extent(MPI_INT, &lb, &(extent[0]));
    MPI_Type_get_extent(MPI_DOUBLE, &lb, &(extent[1]));
    
    offsets[0] = 0;
    offsets[1] = 3 * extent[0];
    offsets[2] = offsets[1] + 12 * extent[1];
    
    MPI_Type_create_struct(3, blockcounts, offsets, oldtypes, MPI_TNODETYPE);
    MPI_Type_commit(MPI_TNODETYPE);
    
    return;
}





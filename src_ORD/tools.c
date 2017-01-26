/* tool functions for use by treecode routines */
#include <math.h>
#include "array.h"
#include "structs.h"


double minval(double *x, int numels) {
    int i;
    double min;
    min = x[0];
    
    for (i = 1; i < numels; i++)
        if (min > x[i]) min = x[i];
    
    return min;
}



double maxval(double *x, int numels) {
    int i;
    double max;
    max = x[0];
    
    for (i = 1; i < numels; i++)
        if (max < x[i]) max = x[i];
    
    return max;
}



double sum(double *x, int numels) {
    int i;
    double sum = 0.0;
    
    for (i = 0; i < numels; i++)  sum = sum + x[i];
    
    return sum;
}



double max3(double a, double b, double c) {
    double max;
    
    max = a;
    if (max < b) max = b;
    if (max < c) max = c;
    
    return max;
}



double min3(double a, double b, double c) {
    double min;
    
    min = a;
    if (min > b) min = b;
    if (min > c) min = c;
    
    return min;
}



int mod(int a, int b)
{
    return (((a) % (b)) + (b)) % (b);
}



void inttobinstr(int value, int length, char* output)
{
    int i;
    output[length] = '\0';
    for (i = length - 1; i >= 0; --i, value >>= 1)
        output[i] = (value & 1) + '0';
}



void remove_particle(struct par *array, int index, int array_length)
{
    int i;
    for(i = index; i < array_length - 1; i++)
        array[i] = array[i + 1];
    
}



void minmaxpar(struct par *parlist, int numpars, double *pminmax) {
    int i;
    
    for (i = 0; i < 3; i++) pminmax[2*i]   = parlist[0].r[i];
    for (i = 0; i < 3; i++) pminmax[2*i+1] = parlist[0].r[i];
    
    for (i = 1; i < numpars; i++) {
        pminmax[0] = fmin(pminmax[0], parlist[i].r[0]);
        pminmax[1] = fmax(pminmax[1], parlist[i].r[0]);
        pminmax[2] = fmin(pminmax[2], parlist[i].r[1]);
        pminmax[3] = fmax(pminmax[3], parlist[i].r[1]);
        pminmax[4] = fmin(pminmax[4], parlist[i].r[2]);
        pminmax[5] = fmax(pminmax[5], parlist[i].r[2]);
    }
    
    return;
}


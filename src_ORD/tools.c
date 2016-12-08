/* tool functions for use by treecode routines */
#include <math.h>
#include <stdio.h>
#include "array.h"
#include "structs.h"


double minval(double *x, int numels)
{
        double min;

        min = x[0];

        for (int i = 1; i < numels; i++) {
                if (min > x[i])
                        min = x[i];
        }

        return min;
}




double maxval(double *x, int numels)
{
        double max;

        max = x[0];

        for (int i = 1; i < numels; i++) {
                if (max < x[i])
                        max = x[i];
        }

        return max;
}




double sum(double *x, int numels)
{
        double sum = 0.0;

        for (int i = 0; i < numels; i++) 
                sum = sum + x[i];

        return sum;
}




double max3(double a, double b, double c)
{
        double max;

        max = a;

        if (max < b)
                max = b;

        if (max < c)
                max = c;

        return max;
}




double min3(double a, double b, double c)
{
        double min;

        min = a;

        if (min > b)
                min = b;

        if (min > c)
                min = c;

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


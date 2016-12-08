/* tool functions for use by treecode routines */
#ifndef H_TOOLS_H                                                                     
#define H_TOOLS_H

#include "structs.h"

double minval(double *x, int numels);
double maxval(double *x, int numels);

double sum(double *x, int numels);

double max3(double a, double b, double c);
double min3(double a, double b, double c);

int mod(int a, int b);

void inttobinstr(int val, int len, char* output);

void remove_particle(struct par *array, int index, int array_length);

#endif /* H_TOOLS_H */

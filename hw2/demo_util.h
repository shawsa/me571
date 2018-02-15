#ifndef DEMO_UTIL_H
#define DEMO_UTIL_H

#include <time.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

/* Array handling routines.  These allocate or delete memory */
void empty_array(int n,double **x);
void ones_array(int n,double **x);
void random_array(int n, double **array);
void linspace_array(double a,double b,int n,double **x);
void delete_array(double **x);

/* Operations on arrays */
double sum_array(int n, double *x);

/* I/O routines */
void read_int(int argc, char** argv, char arg[], int* value, int *err);
void print_global(const char* format, ... );
void print_debug(const char* format, ... );

/* Miscellaneous */
void set_rank();
double random_number();
void random_seed();
int pow2(int p);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif

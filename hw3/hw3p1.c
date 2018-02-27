# include <math.h>
# include <mpi.h>
# include "demo_util.h"
# include <stdio.h>



int main(int argc, char** argv){  

    int rank, nprocs;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* ----------------------------------------------------------------
       Read parameters from the command line
     ---------------------------------------------------------------- */
    int p, itermax;
    double tol;
    int err;
    if (rank == 0)
    {        
        int m,err,loglevel;
        read_int(argc,argv, "-p", &p, &err);        
        read_int(argc,argv, "--itermax", &itermax, &err);
        read_double(argc,argv, "--tol", &tol, &err);   

        printf("p: %d\n", p);
        printf("itermax: %d\n", itermax);
        printf("tol: %g\n", tol); 
    }
    
    MPI_Finalize();
    return 0;
}

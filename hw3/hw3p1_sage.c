# include <math.h>
# include <mpi.h>
# include "demo_util.h"
# include <stdio.h>

/********************************************************************

    f(x) denotes the forcing term
    Au=b is the linear system we're solving
    M = 2^p is the number of unknowns per process
    N = M*procs is the total number of unkowns (or interior nodes)

    Each process will have M+2 nodes stored. M interior nodes and 
        2 end point nodes that it will synchronize at each iteration.

    Residuals and updates will be computed on the fly to conserve memory.

    The convention will be that each process recieves their left endpoint
        from the previous process before passing on their right end point
        to the next process.

********************************************************************/

double foo(double x){
    return -1*(2*PI)*(2*PI)*cos(2*PI*x);
}


int main(int argc, char** argv){  

    int rank, nprocs;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Status status;

    int p, M, itermax;
    double tol;
    int err;
    if (rank == 0)
    {        
        int m,err,loglevel;
        read_int(argc,argv, "-p", &p, &err);        
        read_int(argc,argv, "--itermax", &itermax, &err);
        read_double(argc,argv, "--tol", &tol, &err);   

        M = pow2(p);

        printf("p: %d\n", p);
        printf("itermax: %d\n", itermax);
        printf("tol: %g\n", tol); 
    }
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD); 
    printf("my M: %d\n", M);

    //Determine interval size
    double subint_size = 1/nprocs;
    double a = subint_size*rank;
    double b = a + subint_size;
    double h = (b-a)/(M+1);

    //Initialize empty vectors

    double u[M+2];
    int i;
    for(i=0; i<M+2; i++){u[i]=0;}

    //Assign end conditions
    if(rank==0){
        u[0] = 1;
    }else if(rank==nprocs-1){
        u[M+1] = 1;
    }

    //Iterate
    int iter;
    double diff, largest_diff, ri, xi;
    for(iter=0; iter<itermax; iter++){
        //hold replaced value
        u_minus_old = u[0];
        xi = a
        //perform a Jacobi iteration, keeping track of largest update difference
        for(i=1; i<M+1; i++){
            xi += h;
            ri = u_minus_old - 2*u[i] + u[i+1] - h*h*foo(xi)
            u_minus_old = u[i];
            u[i] += 0.5*ri;
            diff = abs(u[i] - u_minus_old);
            if(diff>largest_diff){largest_diff = diff;}
        }
        //find largest update difference and break if below tol
        if(largest_diff<tol){break;}

        //Synchronize end conditions
        if(rank!=0){
            //receive left end condition
        }
        if(rank != nprocs-1){
            //send right end condition
            
        }
    }

    //print output in order
    
    MPI_Finalize();
    return 0;
}

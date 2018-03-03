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
    return -1*(2*M_PI)*(2*M_PI)*cos(2*M_PI*x);
}


int main(int argc, char** argv){  

    int rank, nprocs;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Status status;

    int p, M, itermax;
    itermax = 10;
    double tol;
    int err;
    if (rank == 0)
    {        
        int m,err,loglevel;
        read_int(argc,argv, "-p", &p, &err);        
        read_int(argc,argv, "--itermax", &itermax, &err);
        read_double(argc,argv, "--tol", &tol, &err);   

        M = pow2(p);

    }
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD); 
    MPI_Bcast(&itermax, 1, MPI_INT, 0, MPI_COMM_WORLD);
    //printf("itermax %d\n", itermax);
    MPI_Bcast(&tol, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //printf("my M: %d\n", M);
    //printf("itermax: %d\n", itermax);
    //printf("tol: %g\n", tol);

    //Determine interval size
    double subint_size = 1.0/nprocs;
    double a = subint_size*rank;
    double b = a + subint_size;
    double h = (b-a)/(M+1);

    //Initialize empty vectors

    double u[M+2];
    int i;
    for(i=0; i<M+2; i++){u[i]=-.3;}

    //Assign end conditions
    if(rank==0){
        u[0] = 1;
    }else if(rank==nprocs-1){
        u[M+1] = 1;
    }

    //Iterate
    int iter;
    double u_minus_old, diff, largest_diff, ri, xi;
    for(iter=0; iter<itermax; iter++){
        //printf("test\n");
        //hold replaced value
        u_minus_old = u[0];
        xi = a;
        //perform a Jacobi iteration, keeping track of largest update difference
        for(i=1; i<M+1; i++){
            xi += h;
            ri = u_minus_old - 2*u[i] + u[i+1] - h*h*foo(xi);
            u_minus_old = u[i];
            u[i] += 0.5*ri;
            diff = abs(u[i] - u_minus_old);
            if(diff>largest_diff){largest_diff = diff;}
        }
        //find largest update difference and break if below tol
        //if(largest_diff<tol){break;}

        //Synchronize end conditions
        if(rank!=0){
            //receive left boundary
            MPI_Recv(&u[0], 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //send left boundary
            MPI_Send(&u[1], 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
        }
        if(rank != nprocs-1){
            //send right boundary
            MPI_Send(&u[M], 1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
            //recieve left boundary
            MPI_Recv(&u[M+1], 1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    //print output in order
        if(rank==0){printf("%g\n",u[0]);}
        if(rank!=0){
            //wait for signal
            MPI_Recv(&u[0], 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        //print array
        for(i=1; i<M+1; i++){
            printf("%g\n",u[i]);
        }
        if(rank != nprocs-1){
            
            //send signal
            MPI_Send(&u[M+1], 1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
        }
        if(rank==nprocs-1){printf("%g\n",u[M+1]);}
    //printf("output test: %g\n", u[1]);
    
    MPI_Finalize();
    return 0;
}

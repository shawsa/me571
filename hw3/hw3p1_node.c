# include <math.h>
# include <mpi.h>
# include "demo_util.h"
# include <stdio.h>

/********************************************************************

    Donna's Node Centered Approach

    f(x) denotes the forcing term
    Au=b is the linear system we're solving
    M = 2^p is the number of unknowns per process
    N = M*procs is the total number of unkowns (or interior nodes)

    Each process will have M+2 nodes stored. M interior nodes and 
        2 end point nodes that it will synchronize at each iteration.
        Except the last process. The last process will have M-2 interior
        nodes making the total number of distinct nodes (including the
        boundary) M*nprocs - a power of 2.

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
    //MPI_Status status;

    int p, M, itermax;
    itermax = 10;
    double tol;
    int err;
    if (rank == 0)
    {        
        read_int(argc,argv, "-p", &p, &err);        
        read_int(argc,argv, "--itermax", &itermax, &err);
        read_double(argc,argv, "--tol", &tol, &err);   

        M = pow2(p);

    }
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD); 
    MPI_Bcast(&itermax, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tol, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    

    //Determine interval size
    double h = 1.0/(M*nprocs);
    double subint_size = (M)*h;
    double a = M*h*rank;    

    //Initialize empty vectors
    
    //last proc has one fewer nodes
    if(rank==nprocs-1){M-=1;}

    double u[M+2];
    int i;
    for(i=0; i<M+2; i++){u[i]=0;}

    //Assign end conditions
    if(rank==0){
        u[0] = 1;
    }
    if(rank==nprocs-1){
        u[M+1] = 1;
    }

    //Iterate
    int iter;
    double u_minus_old, diff, largest_diff_global, largest_diff, ri, xi;
    for(iter=0; iter<itermax; iter++){
        largest_diff = 0;
        //hold replaced value
        u_minus_old = u[0];
        xi = a;
        //perform a Jacobi iteration, keeping track of largest update difference
        for(i=1; i<M+1; i++){
            xi += h;
            ri = u_minus_old - 2*u[i] + u[i+1] - h*h*foo(xi);
            u_minus_old = u[i];
            u[i] += 0.5*ri;
            diff = fabs(u[i] - u_minus_old);
            if(diff>largest_diff){largest_diff = diff;}
        }

        //find largest update difference and break if below tol
        MPI_Reduce(&largest_diff, &largest_diff_global, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Bcast(&largest_diff_global, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
        if(largest_diff_global < tol){break;}

        //Synchronize end conditions
        if(rank!=0){
            //receive left boundary
            MPI_Recv(&(u[0]), 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //send left boundary
            MPI_Send(&(u[1]), 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
        }
        if(rank != nprocs-1){
            //send right boundary
            MPI_Send(&(u[M]), 1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
            //recieve right boundary
            MPI_Recv(&u[M+1], 1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    //print output in order
        if(rank==0){
            printf("%.19g\n",u[0]);
        }
        if(rank!=0){
            //wait for signal
            MPI_Recv(&u[0], 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        //print array
        for(i=1; i<M+1; i++){
            printf("%.19g\n",u[i]);
        }
        if(rank != nprocs-1){
            
            //send signal
            MPI_Send(&(u[M+1]), 1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
        }
        if(rank==nprocs-1){
            printf("%.19g\n",u[M+1]);
            printf("%d\n", iter);
            printf("%.19g\n",largest_diff_global);
        }
    //printf("output test: %g\n", u[1]);
    
    MPI_Finalize();
    return 0;

    /*
        printf("rank: %d\n",rank);
        for(i=0; i<M+1; i++){
            printf("%g\t", u[i]);
        }
        printf("\n");
    */
    
}

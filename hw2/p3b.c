# include <math.h>
# include <mpi.h>
# include "demo_util.h"
# include <stdio.h>


double foo(double x){
    return -1*(2*M_PI)*(2*M_PI)*sin(2*M_PI*x);
}

double exact(double x){
    return sin(2*M_PI*x);
}

double goo(double x, double t){
    if(t<x){
        return (x-1)*t*foo(t);
    }else{
        return (t-1)*x*foo(t);
    }
}

double trap_rule(double t0, double t1, double x, long n){
    double h = (t1 - t0)/n;
    double ret = 0.5*(goo(x, t0)+goo(x,t1));
    double step = (t1-t0)/n;
    int i;
    for(i=1; i<n; i++){
        ret += goo(x,t0+step*i);
    }
    return ret*h;
}


int main(int argc, char** argv){

    double a = 0;
    double b = 1;
    double trap_total;
    double* u;    

    int rank, nprocs;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    //read arguments and broadcast
    long N, n;
    if(rank==0){
        int p, err;
        read_int(argc,argv, "-p", &p, &err);
        N = pow(2,p);
    }
    MPI_Status status;
    MPI_Bcast(&N, 1, MPI_LONG, 0, MPI_COMM_WORLD);

    n = N/nprocs;

    int i;
    double xi;
    double h = (b-a)/N;
    for(i=0; i<=N; i++){
        xi = a + h*i;
        double sub_int = (b-a)/nprocs;
        double t0 = a + sub_int*rank;
        double t1 = t0 + sub_int;
        double trap = trap_rule(t0, t1, xi, n);

        MPI_Reduce(&trap, &trap_total, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if(rank==0){
            printf("%.19g, %.19g, %.19g\n", xi, trap_total, trap_total-exact(xi));
        }
    }
    
    MPI_Finalize();
    return 0;
}

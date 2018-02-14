# include <math.h>
# include <mpi.h>
# include "demo_util.h"
# include <stdio.h>


double foo(double x){
    return (x-1)*(x-1)*exp(-1*x*x);
}

double trap_rule(double x0, double x1, int n){
    double h = (x1 - x0)/n;
    double ret = 0.5*(foo(x0)+foo(x1));
    double step = (x1-x0)/n;
    for(int i=1; i<n; i++){
        ret += foo(x0+step*i);
    }
    return ret*h;
}


int main(int argc, char** argv){

    double a = -1;
    double b = 1;
    double trap_total;
    
    int rank, nprocs;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    //read arguments and broadcast
    int n;
    if(rank==0){
        int p,N, err;
        read_int(argc,argv, "-p", &p, &err);
        printf("p = %d\n", p);
        N = pow(2,p);
        printf("N = %d\n", N);
        n = N/nprocs;
    }
    MPI_Status status;
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);    

    double sub_int = (b-a)/nprocs;
    double x0 = a + sub_int*rank;
    double x1 = x0 + sub_int;
    double trap = trap_rule(x0, x1, n);

    MPI_Reduce(&trap, &trap_total, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(rank==0){
        double exact = 1.872592957265838754602878538234098148617687929406051152575;
        printf("%g\n", trap_total - exact);
    }
    MPI_Finalize();
    return 0;
}

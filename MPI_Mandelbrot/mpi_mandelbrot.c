#include <mpi.h>
#include <stdio.h>


char gen_pixel(double r, double i, int thresh){
	double zr = 0;
	double zi = 0;
	double new_zr;
	double new_zi;
	int j = 0;
	for(j=0; j<thresh; j++){
		new_zr = zr*zr - zi*zi + r;
		new_zi = 2*zr*zi + i;
		zr = new_zr;
		zi = new_zi;
		if(zi*zi + zr*zr > 4){break;}
	}
	if(j==thresh){return '@';}
	else{return ' ';}
}

void print_img(int rows, int cols, char img[rows][cols]){
	for(int r=0; r<rows; r++){
		for(int c=0; c<cols; c++){
			printf("%c", img[r][c]);
		}
		printf("\n");
	}
}


void main(int argc, char** argv){

	char img[200][200];
	
	MPI_Init(&argc, &argv);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	printf("Rank %d starting...\n", rank);
	int start = rank*50;
	int end = start + 50;
	for(int row=start; row<end; row++){
		for(int col=0; col<200; col++){
			double x = -2 + col*(4.0/200.0);
			double y = 2 - row*(4.0/200.0);
			img[row][col] = gen_pixel(x,y,50);
		}
	}

	MPI_Finalize();
	print_img(200, 200, img);
}

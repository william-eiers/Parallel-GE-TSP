#include <iostream>
#include "Chromosome.h"
#include "World.h"
#include "mpi.h"
#include "omp.h"


using namespace std;

// to run program...
//     ./<executable> <initial population> <number of generations> <distance file>
//     				  <mutation chance> <parent/child ratio> <number of cities>
//     				  <read/random option>
int main(int argc, char** argv) {
	
     
	MPI_Init(&argc, &argv);
	int myrank, nprocs;
        MPI_Comm_rank(MPI_COMM_WORLD,&myrank) ;
        MPI_Comm_size(MPI_COMM_WORLD,&nprocs) ;
        srand((myrank+1)*time(NULL)) ;
        
	World world(atoi(argv[1]), atoi(argv[2]), argv[3], atoi(argv[4]), atoi(argv[5]),
				atoi(argv[6]), atoi(argv[7]));
	world.simulate();
     
	// gather the best result from everyone, output the best among those
        
	if(myrank == 0) {
		int numcities = world.getNumberCities();
		// holding the best tour
		int *bestTour = new int[numcities];
		memcpy(bestTour, world.findFittest()->getTour(), numcities*sizeof(int));
		double bestDist = world.findFittest()->getFitness();
		// buffer for receiving the best tour
		int *buffer = new int[numcities];
		double dist;
		for(int i = 1; i < nprocs; i++) {
			// receive tour and distance from every processor
			MPI_Recv(buffer, numcities, MPI_INT, i,i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			MPI_Recv(&dist, 1, MPI_DOUBLE, i, i*nprocs, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			// if tour received is better, make it the best one
                        
			if(dist < bestDist) {
				memcpy(bestTour,buffer,numcities*sizeof(int));
				bestDist = dist;
			}
		}
		cout << "Best Tour: ";
		for(int i = 0; i < numcities; i++) {
			cout << bestTour[i] << " ";
		}
		cout << endl << "Distance: " << bestDist << endl;
        cout << "MPI used: " << nprocs << " procs" << endl;
        cout << "OpenMP used " << world.numProcs << " procs" << endl;
        delete[] bestTour;
        delete[] buffer;
	// every other processer sends best tour to processor 0
	} else {
		double dist = world.findFittest()->getFitness();
                
		MPI_Send(world.findFittest()->getTour(), world.getNumberCities(), MPI_INT, 0, myrank, MPI_COMM_WORLD);
		MPI_Send(&dist, 1, MPI_DOUBLE, 0, myrank*nprocs, MPI_COMM_WORLD);
                
	}
         
	MPI_Finalize();
	
	return 0;
}

#include "World.h"

World::World(unsigned num_pop, unsigned max_generations, char* distancefile,
			 unsigned mutation_chance, unsigned child_ratio, unsigned x, unsigned option) {
	this->num_pop = num_pop;
	this->max_generations = max_generations;
	this->mutation_chance = mutation_chance;
	this->child_ratio = child_ratio;
	// if option 1, read from file
	if(option == 1) {
		population = new vector<Chromosome*>();
		num_cities = read(distancefile);
	// if option 2, generate random data
	} else {
		population = new vector<Chromosome*>();
		num_cities = randomizeCities(x,1100,4000);
	}
}

World::~World() {
	for(int i = 0; i<population->size(); i++) {
		delete (*population)[i];
	}
	delete population;
	
	for(int i = 0; i < num_cities; i++) {
		delete[] dmatrix[i];
	}
	delete[] dmatrix;
}

// continually evolve the population from generation
// to generation, until a "good enough" solution is found
// or until max generations has been reached.
void World::simulate() {
    int myrank;
    int nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    // how much to slice off of the current population each generation
	double sliceRatio = 0.00;
	// how many children per generation get produced
	int nextPop,max_child;
    // mutation chance, should be between 5-10%
	double mChance = mutation_chance * (myrank+1);
    double maxSimilarity; 0.09;
    int initialPop = num_pop;
	int max_parents;
	// create the initial random population
	for(unsigned i = 0; i < num_pop; i++) {
		population->push_back(new Chromosome(num_cities, dmatrix));
	}
	int threads[8] = {0};
	// evolve current population, for max_generations generations
    numProcs = omp_get_num_procs();
    for(unsigned l = 0; l < max_generations; l++) {
        nextPop = num_pop - num_pop*sliceRatio;
        maxSimilarity = maxSimilarity - maxSimilarity*sliceRatio;
        max_child = nextPop / child_ratio;
		max_parents = nextPop - max_child;
		vector<Chromosome*> *children = new vector<Chromosome*>(max_child,NULL);
		// use tournament selection to select two random parents
		// use random parents to create new child
		// add child to child population
                
        // use shared-memory OpenMP to speed up making children
        #pragma omp parallel
        {
			#pragma omp for
            for(int i = 0; i < max_child; i++) {

                const Chromosome *p1 = select(5);
                const Chromosome *p2 = select(5);
                // TODO : check duplicates
                double chance1 = p1->similarity(*p2);
                //Review this with Will. Results are appreantly better
                (*children)[i] = (chance1 < maxSimilarity) ? new Chromosome(*p1, *p2) : new Chromosome(*(select(5)), *(select(5)));
                 // with % chance to mutate...
                if(((double) rand() / (RAND_MAX)) <= mChance)
                    (*children)[i]->mutate();
            }
        }
                // use gnu parallel to sort the population in parallel
		__gnu_parallel::sort(population->begin(), population->end(), compareChromos);
        // take a good portion of the best parents
		
        
        int max_best = (int)(max_parents * .6); // SUBJECT TO CHANGE
        // take a medium portion of the mediocre parents
        int max_medium = (int)(max_parents * .3); // SUBJECT TO CHANGE
        // take a low portion of the worst parents
		int max_worst = (int)(max_parents * .1); // SUBJECT TO CHANGE
		int rand_index;
		// select some of the parents to continue on
		for(int i = 0; i < max_best; i++) {
			rand_index = (rand() % (int)(0.33 * num_pop));
			children->push_back(new Chromosome(*((*population)[rand_index])));
		}
        
		for(int i = 0; i < max_medium; i++) {
			rand_index = (int)(rand() % (int)(0.33 * num_pop) + 0.33*num_pop);
			children->push_back(new Chromosome(*((*population)[rand_index])));
		}
		for(int i = 0; i < max_worst; i++) {
			rand_index = (int)(rand() % int(0.33 * num_pop) + 0.66*num_pop);
			children->push_back(new Chromosome(*((*population)[rand_index])));
		}
		
		// delete any allocated memory
		for(unsigned i = 0; i < population->size(); i++) {
			delete (*population)[i];
			population->pop_back();
		}
		// clear previous population
		// set new population as the children + parents from last generation
		num_pop = children->size();
		population->swap(*children);
		delete children;
		children = NULL;
		migrate(myrank, nprocs);
	}
}

void World::migrate(int myrank, int nprocs) {
    if(nprocs == 1) return;
    best = findFittest();
	// each processor broadcasts their best tour
    int *tours = new int[nprocs*num_cities];
    double *dbuffers = new double[nprocs];
    
    // use MPI broadcast to give out the tours and distances
    dbuffers[myrank] = best->getFitness();
    
    memcpy(&tours[myrank*num_cities],best->getTour(), num_cities*sizeof(int));
    
    MPI_Allgather(&tours[myrank*num_cities],num_cities, MPI_INT, tours,num_cities, MPI_INT, MPI_COMM_WORLD) ;
    MPI_Allgather(&(dbuffers[myrank]), 1, MPI_DOUBLE, dbuffers, 1, MPI_DOUBLE, MPI_COMM_WORLD) ;
    
    
    // delete the last nprocs-1 chromosomes from array
    for(int i = 0; i < nprocs-1; i++) {
        delete (*population).back();
        population->pop_back();
    }
    
    // add the new tours received
    for(int i = 0; i < nprocs; i++) {
        if(i == myrank) continue;
        population->push_back(new Chromosome(&tours[i*num_cities],num_cities,dbuffers[i],dmatrix));
    }
    
    // delete buffers
    
    delete[] tours;
    delete[] dbuffers;
}

// k-tournament selection
const Chromosome * World::select(int number) const {
	// best chromosome
	Chromosome *b = NULL;
	// how many in tournament
	//int k = 2; // subject to change based on results
	for(int i = 0; i < number; i++) {
		int r = rand() % num_pop;
		if(b == NULL || b->getFitness() > (*population)[r]->getFitness()) {
			b = (*population)[r];
		}
	}
	return b;
}

bool World::compareChromos(Chromosome *a, Chromosome *b) {
	return a->getFitness() < b->getFitness() ? true : false;
}

 int World::read(char *filename) {
	 ifstream inFile;     //Input file stream object
	 string line;
	 int size;
	 int row=0;
	 //Open the file
	 inFile.open(filename);

	 if (inFile.fail())
	 {
		 cout << "Error opening data file!\n";
		 exit(102);
	 }

	 while (inFile.good()){
		 getline(inFile, line);
		 istringstream myStream(line);
		 myStream >> size;
		 //Allocating the 2d array
		dmatrix = new double*[size];
		 for (int i = 0; i < size; ++i)
			dmatrix[i] = new double[size];

		 while (getline(inFile, line)){
			 istringstream myStream(line);
			 int col = 0;
			 while (myStream >> dmatrix[row][col++]);
			 row++;
		 }


	}
	return size;
 }

Chromosome* World::findFittest(){
	 if (population->empty())
		 return NULL;
	 else
	 {
		__gnu_parallel::sort(population->begin(), population->end(), compareChromos);
		 return (*population)[0];
	 }

}

int World::randomizeCities(int x, double fMin, double fMax){

	 dmatrix = new double*[x];
	 for (int i = 0; i < x; ++i)
	 dmatrix[i] = new double[x];

	 for (int i = 0; i < x; ++i)
		 dmatrix[i][i] = 0;
	 for (int i = 0; i < x; ++i)
		 for (int j = i + 1; j < x; j++){
		 double f = (double)rand() / RAND_MAX;
		  dmatrix[i][j] =  fMin + f * (fMax - fMin);
		  dmatrix[j][i] =  fMin + f * (fMax - fMin);
		 }

	 return x;


}

// convert contents to buffer for MPI sending
void World::convertBuffer(const Chromosome* chromo, int *buffer, double &dBuff) const {
	dBuff = chromo->getFitness();
	buffer = new int[num_cities];
	memcpy(buffer, chromo->getTour(), num_cities*sizeof(int));	
} 

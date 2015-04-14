/*
 * CS140 Final Project
 * William Eiers 7859812
 * Armin Akhbari 8285298
 */

#ifndef WORLD_H
#define WORLD_H
//#include "City.h"

#include <sys/time.h>
#include "omp.h"
#include <parallel/algorithm>
#include "mpi.h"
#include <string.h>
#include "Chromosome.h"
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>

using namespace std;

class World {
	vector<Chromosome*> *population;
	Chromosome *best;
	unsigned num_pop;
	unsigned max_generations;
	unsigned num_cities;
	unsigned mutation_chance;
	unsigned child_ratio;
	double **dmatrix;
        
        
        
	// tournament selection
	const Chromosome * select(int number) const;
	static  bool compareChromos(Chromosome *a, Chromosome *b);

	//input data
	int read(char *);
	int randomizeCities(int, double, double);
	// private, useless
	World() {}

public:
	// constructor for num_pop, max_gen, distancefile, mutation_chance and child_ratio
	// number of cities, option (read/random)
	World(unsigned, unsigned, char*, unsigned, unsigned, unsigned, unsigned);
	~World();

	// get cities from database
	void create();
	// generate distance matrix based upon cities
	void genMatrix();
	// workhorse of the algorithm
	void simulate();
	// migrate the best between processors
	void migrate(int,int);
	Chromosome* findFittest();
	const Chromosome* getBest() const { return (best); }
	void convertBuffer(const Chromosome* , int *, double &) const;

	//getters and setters
	void setNumberPopulation(int num) { num_pop = num; }
	int getNumberPopulation() const { return num_pop; }
	void setMaxGenerations(int gen) { max_generations = gen; }
	int getMaxGenerations() const { return max_generations; }
	void setNumberCities(int num) { num_cities = num; }
	int getNumberCities() const { return num_cities; }

	int numProcs;
};


#endif

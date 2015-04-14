/*
 * CS140 Final Project
 * William Eiers 7859812
 * Armin Akhbari 8285298
 */

#ifndef CHROMOSOME_H
#define CHROMOSOME_H
//#include "City.h"

#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
using namespace std;
class Chromosome{
	int* tour; 	                                     //dna representation
	unsigned N;											 //number of cities
	double tourLength;                               //fitness score
	double **dptr;									 //pointer to distance matrix (N*N)
	
	void calcFitness();
	int getCityIndex(int) const;
	
public:
	Chromosome(unsigned, double**); 	//check                 //empty chromosome, will fill with random tour
	Chromosome(int*, unsigned, double, double **);              //fill with
	Chromosome(const Chromosome &);
	Chromosome(const Chromosome &, const Chromosome &);          //crossover, create new chromo
	~Chromosome();
	void mutate(); //check
	int* getTour() const {return tour;}		 //check
	
	// getters & setters
	unsigned getCityCount() const {return N;}
	void setFitness(double fitness) {tourLength = fitness;} //check
	double getFitness() const {return tourLength;} //check
	void setDistanceMatrix(double** d) {dptr = d;} //check
	double** getDistanceMatrix() const{return dptr;} //check
	void showTour();
	double similarity(const Chromosome&) const;

};

#endif 

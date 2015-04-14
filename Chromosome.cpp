/*
 * CS140 Final Project
 * William Eiers 7859812
 * Armin Akhbari 8285298
 */

#include "Chromosome.h"
void Chromosome::calcFitness() {

	if (dptr == NULL)
	{
		cout << "Your distance Matrix is NULL\n";
		exit(-1);
	}
	tourLength = 0;
	for (unsigned i = 0; i < N - 1; i++)
		tourLength += dptr[tour[i]][tour[i+1]];
		
	tourLength += dptr[tour[N - 1]][tour[0]];
}

// finds the index of a city in the tour
// TODO: MAYBE HASHTABLE
int Chromosome::getCityIndex(int city) const {
	for(unsigned i = 0; i < N; i++) {
		if(tour[i] == city) return i;
	}
	return -1;
}

//random tour chromosome
Chromosome::Chromosome(unsigned n, double** distances) {
	N = n;
	tourLength = 0;
	dptr = distances;
	tour = new int[N];
	for(unsigned i = 0; i < N; i++) {
		tour[i] = i;
	}
	
	// randomize by swapping random indices
	random_shuffle(&tour[0], &tour[N-1]);

	calcFitness();
}

Chromosome::Chromosome(int *t, unsigned n, double distance, double **distances) {
	N = n;
	dptr = distances;
	tourLength = distance;
	tour = new int[N];
	memcpy(tour, t, N*sizeof(int));
}

Chromosome::Chromosome(const Chromosome &chromo) {
	N = chromo.N;
	dptr = chromo.dptr;
	tourLength = chromo.tourLength;
	tour = new int[N];
	for (unsigned i = 0; i < N; i++) {
		tour[i] = chromo.tour[i];
	}
}

// greedy crossover/recombination to create a new chromo
Chromosome::Chromosome(const Chromosome &ch1, const Chromosome &ch2) {
	// basic setup
	N = ch1.N;
	tour = new int[N];
	dptr = ch1.dptr;
	tourLength = 0;
	
	// recombination
	// select first city from first chromo
	// check if either a duplicate city
	// 		if yes, take other city
	// 		if both, randomly select non-selected city
	//		if neither, compare distances between next city of both
	//			chromosomes
	//		take the shorter path
	int nextStop = 0;
	int city_ch1 = 0;
	int city_ch2 = 0;
	int next_city_ch1 = 1;
	int next_city_ch2 = 1;
	int lastChromoTaken = 0;
	bool *taken = (bool *)calloc(N, sizeof(bool));
	
	tour[nextStop++] = ch1.tour[0];
	taken[ch1.tour[0]] = true;
	lastChromoTaken = 1;
	// populate tour
	while(nextStop != N) {
		// get next city
		if(lastChromoTaken == 1) {
			city_ch2 = ch2.getCityIndex(ch1.tour[city_ch1]);
			next_city_ch2 = (city_ch2 == N-1 ? 0 : city_ch2 + 1);
		} else if(lastChromoTaken == 2) {
			city_ch1 = ch1.getCityIndex(ch2.tour[city_ch2]);
			next_city_ch1 = (city_ch1 == N-1 ? 0 : city_ch1 + 1);
		} else {
			city_ch1 = getCityIndex(tour[nextStop-1]);
			city_ch2 = getCityIndex(tour[nextStop-1]);
			next_city_ch1 = (city_ch1 == N-1 ? 0 : city_ch1 + 1);
			next_city_ch2 = (city_ch2 == N-1 ? 0 : city_ch2 + 1);
		}
		// determine if either city is taken, following rules above
		// both taken
		if(taken[ch1.tour[next_city_ch1]] && taken[ch2.tour[next_city_ch2]]) {
			// take a random city that hasn't been taken yet
			int r = rand() % N;
			int cnt = r;
			while(taken[cnt] && cnt+1 != r) {
				cnt = (cnt == N-1 ? 0 : cnt+1);
			}
			// if can't find any free cities, catastrophic error
			if(taken[cnt]) {
				cout<<"catastrophic error, Chromosome.cpp, constructor\n";
				exit(-1);
			}
			taken[cnt] = true;
			tour[nextStop++] = cnt;
			tourLength += dptr[tour[nextStop-2]][tour[nextStop-1]];
			lastChromoTaken = 3;
		}
		// only next city of second taken
		else if(taken[ch2.tour[next_city_ch2]]) {
			taken[ch1.tour[next_city_ch1]] = true;
			tour[nextStop++] = ch1.tour[next_city_ch1];
			tourLength += dptr[tour[nextStop-2]][tour[nextStop-1]];
			city_ch1 = next_city_ch1;
			next_city_ch1 = (next_city_ch1 == N-1 ? 0 : next_city_ch1 + 1);
			lastChromoTaken = 1;
		}
		// next city of first taken
		else if(taken[ch1.tour[next_city_ch1]]) {
			taken[ch2.tour[next_city_ch2]] = true;
			tour[nextStop++] = ch2.tour[next_city_ch2];
			tourLength += dptr[tour[nextStop - 2]][tour[nextStop - 1]];
			city_ch2 = next_city_ch2;
			next_city_ch2 = (next_city_ch2 == N-1 ? 0 : next_city_ch2 + 1);
			lastChromoTaken = 2;
		}
		// neither next city taken
		else {
			double d1 = dptr[ch1.tour[city_ch1]][ch1.tour[next_city_ch1]];
			double d2 = dptr[ch2.tour[city_ch2]][ch2.tour[next_city_ch2]];
			// take the first chromosomes city
			if(d1 < d2) {
				taken[ch1.tour[next_city_ch1]] = true;
				tour[nextStop++] = ch1.tour[next_city_ch1];
				tourLength += dptr[tour[nextStop-2]][tour[nextStop-1]];
				city_ch1 = next_city_ch1;
				next_city_ch1 = (next_city_ch1 == N-1 ? 0 : next_city_ch1 + 1);
				lastChromoTaken = 1;
			}
			// take the second chromosomes city
			else {
				taken[ch2.tour[next_city_ch2]] = true;
				tour[nextStop++] = ch2.tour[next_city_ch2];
				tourLength += dptr[tour[nextStop - 2]][tour[nextStop - 1]];
				city_ch2 = next_city_ch2;
				next_city_ch2 = (city_ch2 == N-1 ? 0 : city_ch2 + 1);
				lastChromoTaken = 2;
			}
		}
	}
	tourLength += dptr[tour[N-1]][tour[0]];
	free(taken);
}

void Chromosome::showTour(){

	for (unsigned i = 0; i <N; ++i)
		cout << tour[i]+1 << " ";
}

// take care of garbage collection
Chromosome::~Chromosome() {
	delete[] tour;
}
// mutate the chromosome by swapping two random cities
// if the mutation results in a better tour, keep it
// otherwise, revert back
void Chromosome::mutate() {
	
	double oldLength;
		int r1, r2;
	// low chance of being identical
	// if they are, just retry once more
	// not worth it to keep trying... worth it ...
	
	do {
		r1 = rand() % N;
		r2 = rand() % N;

	} while (r1 == r2);
	
	oldLength = tourLength;
	swap(tour[r1], tour[r2]);
	calcFitness();
	// if mutation didn't help... revert back
	if(tourLength > oldLength) {
		tourLength = oldLength;
		swap(tour[r1], tour[r2]);
	}
	// otherwise, keep the new (better) tour
	return;
}

double Chromosome::similarity( const Chromosome &x) const{

	int similar = 0;
	for (unsigned i = 0; i < N; ++i)
		if (tour[i] == x.getTour()[i])
			similar++;
	
	return similar / N;


}

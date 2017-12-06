/***************************************************************************
 *   copyright Pleuni PENNINGS                                 *
 *   pennings@fas.harvard.edu                                 *
***************************************************************************/

// System wide headers:
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sstream> 
#include <iostream> 
#include <fstream>
#include <vector> 
#include <math.h>
#include <string>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
using namespace std;

//Global variables
double cost; //cost of mutant 
double mu; //mut rate
#define NUMINDS 1000 //fixed population size
int NumMutants; 
int NumWildtype; 
int mut_from_wt; 
int wt_from_mut; 
double freq;//frequency of the mutant 
double expectedfreq; // 
double pi; // heterozygosity
double MeanFitness; //mean fitness of the populations (not sure if needed)


//Variables needed in the program
unsigned int seed; 
int output_every_Xgen; // how often to output data
double numgen_inN; //how many N generations to run the program
int numgen; // how many generations to run the program
int start_output; //when to start writing output
unsigned int repeat; //which of nRuns is running
ofstream output; // to write output

void getParameters(){
	cerr << "seed: "; cin >> seed; 
	cerr << "mu: "; cin >> mu; 
	cerr << "cost: "; cin >> cost; 
	cerr << "op_every_Xgen: "; cin >> output_every_Xgen; 
	cerr << "ngen_inN: "; cin >> numgen_inN; 
	cerr << "st_output: "; cin >> start_output; 
	seed*=100;
	numgen = numgen_inN*NUMINDS;
	cerr << numgen_inN<< " "<<numgen << "\n";
}

void popAdapt(gsl_rng *rng){
	//popAdapt has two parts: initialize (to make a monomorphic population and seed the random number generator) and evolve (a loop that loops until the ancestral allele is lost). 
	//initialize
	cout << "t\tpi\tfreq\n"; //del muts
	NumMutants=0; NumWildtype=NUMINDS;
	//evolve	
	for (int t=0; t<(numgen/output_every_Xgen);t++)  {//each loop is a bunch of generations
		//output pi
		freq=double(NumMutants)/NUMINDS; pi = 2*freq*(1-freq);
		if ((t*output_every_Xgen)>start_output) {
			cout  << t*output_every_Xgen << "\t" << pi << "\t" << freq<< "\n";} 
		
		for (int t2 = 0; t2<(output_every_Xgen); t2++){//each loop is 1 gen
			//I am going to change the order. New order: selection & inheritence before mutation
			//inherit
			freq = double(NumMutants)/double(NUMINDS); //I think that this is the frequency I should output!
			MeanFitness = 1-(freq*cost);
			expectedfreq = (freq * (1-cost)) / MeanFitness;
			NumMutants = gsl_ran_binomial(rng, expectedfreq, NUMINDS); 
			NumWildtype = NUMINDS - NumMutants;
			//mutation
			mut_from_wt = gsl_ran_binomial(rng, mu, NumWildtype);
			wt_from_mut = gsl_ran_binomial(rng, mu, NumMutants);
			NumMutants+=(mut_from_wt-wt_from_mut);
			NumWildtype-=(mut_from_wt-wt_from_mut);
		}
	}
}

int main (int argc, char *argv[]){
	getParameters();
	srand(seed*500000);
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);//choose random number generator
//	cerr << "Seed  " <<seed<< "\n"; 
	gsl_rng_set (rng, (unsigned long)repeat+seed+1); 
//	cerr << "mut" << mu<< "\n";
//	cerr << "cost" << cost<< "\n";
	popAdapt(rng);
	return 5;
}


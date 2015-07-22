/***************************************************************************
© Champak Beeravolu Reddy 2015-now

champak.br@gmail.com

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.

***************************************************************************/

/*
 * main.cpp
 *
 *  Created on: 3 Jun 2015
 *      Author: champost
 */


#include <ctime>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <map>
#include <vector>
#include <algorithm>

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

#include "MersenneTwister.h"
#include "main.h"

using namespace std;

//************ EXTERN **************
int brClass, mutClass, foldBrClass, allBrClasses;
double main_theta, main_rho;
double main_div_time;
//**********************************

MTRand rMT;
map<int, string> mutConfig2Str;
map<int, vector<int> > MutConfig;
map<vector<int>, int> intVec2BrConfig;
vector<int> npopVec;
vector<double> dataTable;
vector<double> finalTable;

int ms_argc;
char **ms_argv;
int ntrees;
bool estimate;
long int finalTableSize;


double ranMT() { return(rMT()); }

void calcFinalTable(double **onetreeTable) {
	for (int i = 0; i < finalTableSize; i++) {
	    double jointPoisson = 1.0;

		for (int j = 0; j < brClass; j++)
			jointPoisson *= onetreeTable[j][MutConfig[i][j]];

		finalTable[i] += jointPoisson;
	}
}

int getBrConfigNum(int *brConfVec) {
	vector<int> vec(brConfVec, brConfVec+npopVec.size());
	return intVec2BrConfig[vec];
}

double computeLik() {
	finalTable = vector<double> (finalTableSize, 0);

	// calling ms
	main_ms(ms_argc, ms_argv);

	double loglik = 0.0;
	if (estimate) {
		for (long int i = 0; i < finalTableSize; i++)
			if (finalTable[i] > 0)
				loglik += log(finalTable[i] / ntrees) * dataTable[i];
	}
	else {
		for (long int i = 0; i < finalTableSize; i++)
//			printf("%.5e\n", finalTable[i]/totProbSum);
			printf("%s : %.5e\n", mutConfig2Str[i].c_str(), finalTable[i]/ntrees);
	}

	return loglik;
}


double optimize_wrapper(const gsl_vector *vars, void *obj) {

	main_theta = gsl_vector_get(vars, 0);
	main_rho = gsl_vector_get(vars, 1);
	main_div_time = gsl_vector_get(vars, 2);
	if ((main_theta < 0) or (main_rho < 0))
		return 999999;
	else
		return -computeLik();
}


//	conversion from decimal to base-mutClass
void evalMutConfigs() {
	int quo, rem;
	for (long int i = 0; i < finalTableSize; i++) {
		quo = i;
		rem = 0;
		stringstream stst;
		stst << ")";
		for (int j = 0; j < brClass; j++) {
			if (quo) {
				rem = quo % mutClass;
				quo /= mutClass;
				if (rem == mutClass-1)
					stst << rem-1 << ">";
				else
					stst << rem;
				MutConfig[i].push_back(rem);
			}
			else {
				stst << "0";
				MutConfig[i].push_back(0);
			}

			if (j < brClass - 1)
				stst << ",";
		}
		stst << "(";
		string config = stst.str();
		reverse(config.begin(),config.end());
		mutConfig2Str[i] = config;
		reverse(MutConfig[i].begin(),MutConfig[i].end());
//		cout << config << endl;
	}
}


void readMutConfigs() {
	string line;
	ifstream ifs("marginals_only.txt",ios::in);
	while (getline(ifs,line))
		dataTable.push_back(atof(line.c_str()));
	ifs.close();
}


void readPopSizes(int npops) {
	string line;
	ifstream ifs("popconfig.txt",ios::in);
	getline(ifs,line);
	stringstream stst;
	stst << line;
	for (int i = 0; i < npops; i++) {
		int tmp;
		stst >> tmp;
		npopVec.push_back(tmp);
	}
	ifs.close();
}


//	conversion from decimal to base-(maxPopSize+1)
void evalBranchConfigs() {
	int quo, rem, maxPopSize, totPopSum, count = 0, sumConfig;
	bool skipConfig;
	maxPopSize = totPopSum = npopVec[0];

	for (size_t i = 1; i < npopVec.size(); i++) {
		totPopSum += npopVec[i];
		if (npopVec[i] > maxPopSize)
			maxPopSize = npopVec[i];
	}
	++maxPopSize;

	for (long int i = 1; i <= (long int) pow(maxPopSize,npopVec.size()); i++) {
		quo = i;
		rem = 0;
		stringstream stst;
		stst << ")";
		sumConfig = 0;
		skipConfig = false;
		vector<int> vec;

		for (size_t j = 0; j < npopVec.size(); j++) {
			if (quo) {
				rem = quo % (maxPopSize);
				quo /= (maxPopSize);
				if (rem > npopVec[npopVec.size()-1-j]) {
					skipConfig = true;
					break;
				}
				stst << rem;
				sumConfig += rem;
				vec.push_back(rem);
			}
			else {
				stst << "0";
				vec.push_back(0);
			}

			if (j < npopVec.size() - 1)
				stst << ",";
		}

		if (sumConfig == totPopSum)
			break;

		if (!skipConfig) {
			stst << "(";
			string config = stst.str();
			reverse(config.begin(),config.end());
			reverse(vec.begin(),vec.end());
			intVec2BrConfig[vec] = count;
			++count;
//			printf("%d\t%d\t%s\n", i, count, config.c_str());
//			printf("%d\t%s\n", count, config.c_str());
		}
	}
}


int main(int argc, char* argv[]) {

	time_t likStartTime, likEndTime;
	int nsam = atoi(argv[1]), kmax = atoi(argv[argc-3]), npopSize = atoi(argv[argc-2]);
	char brFold = argv[argc-1][0];
	ms_argc = argc - 4;
	ms_argv = argv;
	ntrees = atoi(argv[2]);

	readPopSizes(npopSize);

	brClass = npopVec[0]+1;
	for (size_t i = 1; i < npopVec.size(); i++)
		brClass *= (npopVec[i]+1);
	brClass -= 2;
	allBrClasses = brClass;
	foldBrClass = 0;

	//	fold the branch classes
	if (brFold == 'f') {
		if (brClass % 2)
			brClass = (brClass+1) / 2;
		else
			brClass = brClass / 2;
		foldBrClass = 1;
	}

	evalBranchConfigs();

	mutClass = kmax+2;
	finalTableSize = (long int) pow(mutClass, brClass);

	evalMutConfigs();

	if (atoi(argv[argc-4])) {

		estimate = true;

		readMutConfigs();

		const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
		gsl_multimin_fminimizer *s = NULL;
		gsl_vector *ss, *x;
		size_t npar = 3;
		size_t iter = 0, i;
		int status;
		double size;

		/* Initial vertex size vector */
		ss = gsl_vector_alloc (npar);

		/* Set all step sizes to .01 */ //Note that it was originally 1
		gsl_vector_set_all (ss, 1);

		/* Starting point */
		x = gsl_vector_alloc (npar);

		gsl_vector_set (x,0,ranMT()*10);
		gsl_vector_set (x,1,ranMT()*10);
		gsl_vector_set (x,2,ranMT()*10);
		printf ("Start coordinates : \n");
		for (i = 0; i < npar; i++)
			printf ("%.6f ", gsl_vector_get (x, i));
		printf ("\n");

		gsl_multimin_function minex_func;
		minex_func.f = optimize_wrapper;
//		minex_func.params = pt;
		minex_func.n = npar;
		s = gsl_multimin_fminimizer_alloc (T, npar);
		gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
		do {
			//cout<<"Now on iteration "<<iter<<endl;
			iter++;
			status = gsl_multimin_fminimizer_iterate(s);
			if (status!=0) { //0 Means it's a success
				printf ("error: %s\n", gsl_strerror (status));
				break;
			}

			size = gsl_multimin_fminimizer_size (s);
			//status = gsl_multimin_test_size (size, 1e-2);
			status = gsl_multimin_test_size (size, 0.00001); //since we want more precision

			if (status == GSL_SUCCESS) {
				printf ("Converged to a maximum at\n");

				//printf ("%5d ", iter);

				for (i = 0; i < npar; i++)
					printf ("%.6f ", gsl_vector_get (s->x, i));
//				printf ("LnL = %.6f size = %.5e\n", -s->fval, size);
				printf ("LnL = %.6f\n\n", -s->fval);
			}
//			else {
//				printf ("%5d ", iter);
//				for (i = 0; i < npar; i++)
//					printf ("%.6f ", gsl_vector_get (s->x, i));
//				printf ("LnL = %.6f\n", -s->fval);
//			}

		} while (status == GSL_CONTINUE && iter < 1000);

		gsl_vector_free(x);
		gsl_vector_free(ss);
		gsl_multimin_fminimizer_free (s);

//		for (double theta = 2.01; theta < 4.0; theta +=0.25) {
//			for (double rho = 4.01; rho < 6.0; rho +=0.25) {
//
//				main_theta = theta;
//				main_rho = rho;
//
////				printf("%.2f\t%.5e\n",theta, computeLik());
//				printf("%.2f\t%.2f\t%.5e\n",theta, rho, computeLik());
//
//			}	//	rho
//		}	//	theta
	}
	else {
		time(&likStartTime);

		estimate = false;
		computeLik();

		time(&likEndTime);
		printf("\n\nTime taken for computation only : %.5f s", float(likEndTime - likStartTime));

	}

	return 0;
}




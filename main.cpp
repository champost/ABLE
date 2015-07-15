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
int brClass, mutClass;
long int finalTableSize;
double main_theta, main_rho, totSum, *finalTable;
//**********************************

MTRand rMT;
map<int, string> int2string;
map<int, vector<int> > MutConfig;

int ms_argc;
char **ms_argv;
double *dataTable;
bool estimate;


double ranMT() { return(rMT()); }

int getMutConfig(int mutConfNum, int brClassNum) { return MutConfig[mutConfNum][brClassNum]; }

double computeLik() {

	totSum = 0.0;
	for (long int i = 0; i < finalTableSize; i++)
		finalTable[i] = 0.0;

	// calling ms
	main_ms(ms_argc, ms_argv);

	double loglik = 0.0;
	if (estimate) {
		for (long int i = 0; i < finalTableSize; i++)
			if (finalTable[i] > 0)
				loglik += log(finalTable[i] / totSum) * dataTable[i];
	}
	else {
		for (long int i = 0; i < finalTableSize; i++)
//			printf("%.5e\n", finalTable[i]/totSum);
			printf("%s : %.5e\n", int2string[i].c_str(), finalTable[i]/totSum);
	}

	return loglik;
}


double optimize_wrapper(const gsl_vector *vars, void *obj) {

	main_theta = gsl_vector_get(vars, 0);
	main_rho = gsl_vector_get(vars, 1);
	if ((main_theta < 0) or (main_rho < 0))
		return 999999;
	else
		return -computeLik();
}


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
		int2string[i] = config;
		reverse(MutConfig[i].begin(),MutConfig[i].end());
//		cout << config << endl;
	}
}


void readMutConfigs() {
	int count = 0;
	string line;
	ifstream ifs("marginals_only.txt",ios::in);
	while (getline(ifs,line)) {
		dataTable[count] = atof(line.c_str());
		++count;
	}
	ifs.close();
}


int main(int argc, char* argv[]) {

//	time_t likStartTime, likEndTime;
	int nsam = atoi(argv[1]), ntrees = atoi(argv[2]), kmax = 3;
	ms_argc = argc - 1;
	ms_argv = argv;

	brClass = nsam-1;	//4-1=3
	mutClass = kmax+2;	//3+2=5
	finalTableSize = (long int) pow(mutClass,brClass);	//5^3=125

	finalTable = (double *) malloc(finalTableSize * sizeof(double));	//5^3=125
	dataTable = (double *) malloc(finalTableSize * sizeof(double));	//5^3=125

	evalMutConfigs();

	if (atoi(argv[argc-1])) {

		estimate = true;

		readMutConfigs();

		const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
		gsl_multimin_fminimizer *s = NULL;
		gsl_vector *ss, *x;
		size_t npar = 2;
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
				printf ("LnL = %.6f\n", -s->fval);
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
//		time(&likStartTime);

		estimate = false;
		computeLik();

//		time(&likEndTime);
//		printf("\n\nTime taken for computation only : %.5f s", float(likEndTime - likStartTime));

	}

	free(finalTable);
	free(dataTable);

	return 0;
}




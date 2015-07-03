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
//#include <gsl/gsl_sort_double.h>
//#include <gsl/gsl_statistics_double.h>
//#include <gsl/gsl_randist.h>

#include "MersenneTwister.h"
#include "ms_new.h"
//#include "main.h"

extern "C" {
int main_ms(int ms_argc, char *ms_argv[], double ***treeTable);
double *** d3matrix(int x, int y, int z);
void freed3matrix(double ***m, int x, int y);
double ** d2matrix(int x, int y);
void freed2matrix(double **m, int x);
int ** d2int_matrix(int x, int y);
void freed2int_matrix(int **m, int x);

double ranMT();
}

using namespace std;
int treeTableX, treeTableY, treeTableZ;
double main_theta, main_rho;
MTRand rMT;

double ranMT() { return(rMT()); }


int main(int argc, char* argv[]) {

//	time_t likStartTime, likEndTime;
	int nsam = atoi(argv[1]), ntrees = atoi(argv[2]), kmax = 3, count;
	treeTableX = ntrees;	//1000
	treeTableY = nsam-1;	//4-1=3
	treeTableZ = kmax+2;	//3+2=5
	long int finalTableSize = (long int) pow(treeTableZ,treeTableY);	//5^3=125

	double ***treeTable = (double ***) malloc(treeTableX * sizeof(double **)),	//1000*3*5
			*finalTable = (double *) malloc(finalTableSize * sizeof(double)),	//5^3=125
			*dataTable = (double *) malloc(finalTableSize * sizeof(double)),	//5^3=125
			*jointPoisson = (double *) malloc(treeTableX * sizeof(double)),		//1000
			loglik, totSum;
	for (int i = 0; i < treeTableX; i++) {
		treeTable[i] = (double **) malloc(treeTableY * sizeof(double *));
		for (int j = 0; j < treeTableY; j++)
			treeTable[i][j] = (double *) malloc(treeTableZ * sizeof(double));
	}


	map<int, string> int2string;
	map<int, vector<int> > MutConfig;

	int quo, rem;
	for (long int i = 0; i < finalTableSize; i++) {
		quo = i;
		rem = 0;
		stringstream stst;
		stst << ")";
		for (int j = 0; j < treeTableY; j++) {
			if (quo) {
				rem = quo % treeTableZ;
				quo /= treeTableZ;
				stst << rem;
				MutConfig[i].push_back(rem);
			}
			else {
				stst << "0";
				MutConfig[i].push_back(0);
			}

			if (j < treeTableY - 1)
				stst << ",";
		}
		stst << "(";
		string config = stst.str();
		reverse(config.begin(),config.end());
		int2string[i] = config;
		reverse(MutConfig[i].begin(),MutConfig[i].end());
//		cout << config << endl;
	}


	if (atoi(argv[argc-1])) {
		count = 0;
		string line;
		ifstream ifs("marginals_only.txt",ios::in);
		while (getline(ifs,line)) {
			dataTable[count] = atof(line.c_str());
			++count;
		}
		ifs.close();

		for (double theta = 2.01; theta < 4.0; theta +=0.25) {
			for (double rho = 4.01; rho < 6.0; rho +=0.25) {

				loglik = 0.0;
				main_theta = theta;
				main_rho = rho;

				// calling ms
				main_ms(argc-1, argv, treeTable);

				totSum = 0.0;
				for (long int i = 0; i < finalTableSize; i++) {
					finalTable[i] = 0.0;

					for (int trees = 0; trees < treeTableX; trees++) {
						jointPoisson[trees] = 1.0;
						for (int j = 0; j < treeTableY; j++)
							jointPoisson[trees] *= treeTable[trees][j][MutConfig[i][j]];

						finalTable[i] += jointPoisson[trees];
					}
					totSum += finalTable[i];
				}

				for (long int i = 0; i < finalTableSize; i++)
					loglik += log(finalTable[i] / totSum) * dataTable[i];

//				printf("%.2f\t%.5e\n",theta, loglik);
				printf("%.2f\t%.2f\t%.5e\n",theta, rho, loglik);

			}	//	rho
		}	//	theta
	}
	else {
//		time(&likStartTime);

		// calling ms
		main_ms(argc-1, argv, treeTable);

		totSum = 0.0;
		for (long int i = 0; i < finalTableSize; i++) {
			finalTable[i] = 0.0;

			for (int trees = 0; trees < treeTableX; trees++) {
				jointPoisson[trees] = 1.0;
				for (int j = 0; j < treeTableY; j++)
					jointPoisson[trees] *= treeTable[trees][j][MutConfig[i][j]];

				finalTable[i] += jointPoisson[trees];
			}
			totSum += finalTable[i];
		}

		for (long int i = 0; i < finalTableSize; i++)
//			printf("%.5e\n", finalTable[i]/totSum);
			printf("%s : %.5e\n", int2string[i].c_str(), finalTable[i]/totSum);

//		printf("\n%.5f\n", totSum);

		//		time(&likEndTime);
//		printf("\n\nTime taken for computation only : %.5f s", float(likEndTime - likStartTime));

	}

	freed3matrix(treeTable, treeTableX, treeTableY);
	free(finalTable);
	free(dataTable);
	free(jointPoisson);

	return 0;
}




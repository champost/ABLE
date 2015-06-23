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
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_randist.h>

#include "MersenneTwister.h"
#include "ms_new.h"
//#include "main.h"

extern "C" {
int main_ms(int ms_argc, char *ms_argv[], double ****treeTable, int **segs);
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
	int nsam = atoi(argv[1]), ntrees = atoi(argv[2]), kmax = 3;
	long int totTrees;
	treeTableX = ntrees;	//1000
	treeTableY = nsam-1;	//4-1=3
	treeTableZ = kmax+2;	//3+2=5
	size_t finalTableSize = (size_t) pow(treeTableZ,treeTableY);	//5^3=125

	double ****treeTable = (double ****) malloc(treeTableX * sizeof(double ***)),	//1000*X*3*5
			*finalTable = (double *) malloc(finalTableSize * sizeof(double)),	//5^3=125
			*dataTable = (double *) malloc(finalTableSize * sizeof(double)),	//5^3=125
			*jointPoisson = (double *) malloc(treeTableX * sizeof(double)),		//1000
			loglik;
	int **segs = (int **) malloc(treeTableX * sizeof(int *));	//1000*(X+2)

	if (atoi(argv[argc-1])) {
		int count = 0;
		string line;
		ifstream ifs("marginals_only.txt",ios::in);
		while (getline(ifs,line)) {
			dataTable[count] = atof(line.c_str());
			++count;
		}
		ifs.close();

		for (double theta = 0.5; theta < 10.0; theta +=0.5) {
			for (double rho = 0.5; rho < 10.0; rho +=0.5) {

				loglik = 0.0;
				main_theta = theta;
				main_rho = rho;

				totTrees = 0;

				main_ms(argc-1, argv, treeTable, segs);

				for (int trees = 0; trees < treeTableX; trees++)
					totTrees += segs[trees][0];

				for (size_t i = 0; i < finalTableSize; i++)
					finalTable[i] = 0.0;

				count = 0;
				for (int i = 0; i < treeTableZ; i++) {
					for (int j = 0; j < treeTableZ; j++) {
						for (int k = 0; k < treeTableZ; k++) {

							for (int trees = 0; trees < treeTableX; trees++) {
								if (segs[trees][0] > 1) {
									jointPoisson[trees] = 0.0;
									for (int blocks = 0; blocks < segs[trees][0]; blocks++)
										jointPoisson[trees] += treeTable[trees][blocks][0][i]*treeTable[trees][blocks][1][j]*treeTable[trees][blocks][2][k];
//										jointPoisson[trees] += double (segs[trees][blocks+2]/segs[trees][1]) * treeTable[trees][blocks][0][i]*treeTable[trees][blocks][1][j]*treeTable[trees][blocks][2][k];
//										jointPoisson[trees] += double (1/segs[trees][1]) * treeTable[trees][blocks][0][i]*treeTable[trees][blocks][1][j]*treeTable[trees][blocks][2][k];
								}
								else
									jointPoisson[trees] = treeTable[trees][0][0][i]*treeTable[trees][0][1][j]*treeTable[trees][0][2][k];
								finalTable[count] += jointPoisson[trees];
							}
							loglik += log(finalTable[count] / totTrees) * dataTable[count];
							++count;
						}
					}
				}

//				printf("%.2f\t%.5e\n",theta, loglik);
				printf("%.2f\t%.2f\t%.5e\n",theta, rho, loglik);

				for (int i = treeTableX-1; i >= 0; i--)
					freed3matrix(treeTable[i], segs[i][0], treeTableY);
				for (int i = treeTableX-1; i >= 0; i--)
					free(segs[i]);
			}	//	rho
		}	//	theta
		free(treeTable);
		free(segs);
	}
	else {
//		time(&likStartTime);

		totTrees = 0;

		main_ms(argc-1, argv, treeTable, segs);

		for (int trees = 0; trees < treeTableX; trees++)
			totTrees += segs[trees][0];

		for (size_t i = 0; i < finalTableSize; i++)
			finalTable[i] = 0.0;

		int count = 0;
		for (int i = 0; i < treeTableZ; i++) {
			for (int j = 0; j < treeTableZ; j++) {
				for (int k = 0; k < treeTableZ; k++) {

					for (int trees = 0; trees < treeTableX; trees++) {
						if (segs[trees][0] > 1) {
							jointPoisson[trees] = 0.0;
							for (int blocks = 0; blocks < segs[trees][0]; blocks++)
								jointPoisson[trees] += treeTable[trees][blocks][0][i]*treeTable[trees][blocks][1][j]*treeTable[trees][blocks][2][k];
//								jointPoisson[trees] += double (segs[trees][blocks+2]/segs[trees][1]) * treeTable[trees][blocks][0][i]*treeTable[trees][blocks][1][j]*treeTable[trees][blocks][2][k];
//								jointPoisson[trees] += double (1/segs[trees][1]) * treeTable[trees][blocks][0][i]*treeTable[trees][blocks][1][j]*treeTable[trees][blocks][2][k];
						}
						else
							jointPoisson[trees] = treeTable[trees][0][0][i]*treeTable[trees][0][1][j]*treeTable[trees][0][2][k];
						finalTable[count] += jointPoisson[trees];
					}
					printf("(%d,%d,%d) : %.5e\n", i, j, k, finalTable[count]/totTrees);
					++count;
				}
			}
		}

//		time(&likEndTime);
//		printf("\n\nTime taken for computation only : %.5f s", float(likEndTime - likStartTime));

		for (int i = treeTableX-1; i >= 0; i--)
			freed3matrix(treeTable[i], segs[i][0], treeTableY);
		free(treeTable);
		freed2int_matrix(segs,treeTableX);
	}

	free(finalTable);
	free(dataTable);
	free(jointPoisson);

	return 0;
}




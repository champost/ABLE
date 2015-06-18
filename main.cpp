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

#include "ms_new.h"
//#include "main.h"

extern "C" int main_ms(int ms_argc, char *ms_argv[], double ***treeTable);
extern "C" double *** d3matrix(int x, int y, int z);
extern "C" void freed3matrix(double ***m, int x, int y);
extern "C" double ** d2matrix(int x, int y);
extern "C" void freed2matrix(double **m, int x);

using namespace std;
int treeTableX, treeTableY, treeTableZ;
double main_theta, main_rho;
int curve;

int main(int argc, char* argv[]) {

	time_t likStartTime, likEndTime;
	int nsam = atoi(argv[1]), ntrees = atoi(argv[2]), kmax = 3;
	treeTableX = ntrees;	//1000
	treeTableY = nsam-1;	//4-1=3
	treeTableZ = kmax+2;	//3+2=5
	size_t finalTableSize = (size_t) pow(treeTableZ,treeTableY);	//5^3=125

	double ***treeTable = d3matrix(treeTableX, treeTableY, treeTableZ),	//1000*3*5
			*finalTable = (double *) malloc(finalTableSize * sizeof(double)),	//5^3=125
			*dataTable = (double *) malloc(finalTableSize * sizeof(double)),	//5^3=125
			*jointPoisson = (double *) malloc(treeTableX * sizeof(double)),		//1000
			loglik;

	curve = 1;
	if (curve) {
		int count = 0;
		string line;
		ifstream ifs("marginals_only.txt",ios::in);
		while (getline(ifs,line)) {
			dataTable[count] = atof(line.c_str());
			++count;
		}
		ifs.close();

		for (double theta = 0.01; theta < 10.0; theta +=0.05) {
			loglik = 0.0;
			main_theta = theta;

			main_ms(argc, argv, treeTable);

			for (size_t i = 0; i < finalTableSize; i++) {
				finalTable[i] = 0.0;
			}

//			stringstream stst;
//			stst << "Theta_" << theta << ".txt";
//			FILE *pf = fopen(stst.str().c_str(),"w");
			int count = 0;
			for (int i = 0; i < treeTableZ; i++) {
				for (int j = 0; j < treeTableZ; j++) {
					for (int k = 0; k < treeTableZ; k++) {
						for (int trees = 0; trees < treeTableX; trees++) {
							jointPoisson[trees] = treeTable[trees][0][i]*treeTable[trees][1][j]*treeTable[trees][2][k];
							finalTable[count] += jointPoisson[trees];
						}
//						fprintf(pf,"(%d,%d,%d) : %.5e\n", i, j, k, finalTable[count]/ntrees);
						++count;
					}
				}
			}
//			fclose(pf);

			for (size_t i = 0; i < finalTableSize; i++) {
				loglik += log(finalTable[i] / ntrees) * dataTable[i];
//				loglik += finalTable[i] / ntrees * dataTable[i];
			}


			printf("%.2f\t%.5e\n",theta, loglik);
//			printf("%.2f\t%.5e\n",theta, log(loglik));
		}
	}
	else {
		time(&likStartTime);

		main_ms(argc, argv, treeTable);

	//	cout << "\n\n";
	//	for (int i = 0; i < treeTableX; i++) {
	//		printf("Tree : %d\n", i+1);
	//		for (int j = 0; j < treeTableY; j++) {
	//			printf("*** %d-ton branches ***\n", j+1);
	//			if (treeTable[i][j][0] == 1)
	//				printf("Total branch length = 0\n\n");
	//			else {
	//				for (int k = 0; k < treeTableZ; k++)
	//					printf("%5.5lf ", treeTable[i][j][k]);
	//				cout << "\n\n";
	//			}
	//		}
	//		cout << "\n\n";
	//	}

		for (size_t i = 0; i < finalTableSize; i++) {
			finalTable[i] = 0.0;
		}

		int count = 0;
		for (int i = 0; i < treeTableZ; i++) {
			for (int j = 0; j < treeTableZ; j++) {
				for (int k = 0; k < treeTableZ; k++) {
					for (int trees = 0; trees < treeTableX; trees++) {
						jointPoisson[trees] = treeTable[trees][0][i]*treeTable[trees][1][j]*treeTable[trees][2][k];
						finalTable[count] += jointPoisson[trees];
					}
					printf("(%d,%d,%d) : %.5e\n", i, j, k, finalTable[count]/ntrees);
					++count;
				}
			}
		}

		time(&likEndTime);

	//	printf("\n\n");
//		count = 0;
//		for (int i = 0; i < treeTableZ; i++) {
//			for (int j = 0; j < treeTableZ; j++) {
//				for (int k = 0; k < treeTableZ; k++) {
//					printf("(%d,%d,%d) : %.5e\n", i, j, k, finalTable[count]/ntrees);
//					++count;
//				}
//			}
//		}

	//	printf("\n\nTime taken for computation only : %.5f s", float(likEndTime - likStartTime));
	}

	freed3matrix(treeTable, treeTableX, treeTableY);
	free(finalTable);
	free(dataTable);
	free(jointPoisson);

	return 0;
}




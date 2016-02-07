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
 * main.h
 *
 *  Created on: 4 Jun 2015
 *      Author: champost
 */

#ifndef MAIN_H_
#define MAIN_H_

using namespace std;

void free_ms_argv();
void profileLik(vector<double> MLEparVec, double maxLnL);
vector<int> getMutConfigVec(unsigned long int i);
string getMutConfigStr(unsigned long int i);
string getMutConfigStr(vector<int> configVec);
void freePoissonProbs();
void calcFinalTable();
double computeLik();
double optimize_wrapper_nlopt(const vector<double> &vars, vector<double> &grad, void *data);
void readDataConfigs();
void evalBranchConfigs();
void readConfigFile(char* argv[]);


extern "C" {
int main_ms(int ms_argc, char *ms_argv[]);
//double *** d3matrix(int x, int y, int z);
//void freed3matrix(double ***m, int x, int y);
//double ** d2matrix(int x, int y);
void freed2matrix(double **m, int x);
//int ** d2int_matrix(int x, int y);
//void freed2int_matrix(int **m, int x);

double ranMT();
void storePoissonProbs(double **onetreeTable);
int getPopSampleStatus(int pop);
int getBrConfigNum(int *brConfVec);
}

#endif /* MAIN_H_ */

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

const string datestring=__DATE__;
const string timestring=__TIME__;

void free_objects();
void profileLik(vector<double> MLEparVec);
vector<int> getMutConfigVec(unsigned long int i);
string getMutConfigStr(unsigned long int i);
string getMutConfigStr(vector<int> configVec);
int getBrConfigNum(vector<int> brConfVec);
void process_tree_cond_bSFS (double ***onetreePoisTable);
void process_tree_exact_bSFS (double ***onetreePoisTable);
void calcBSFSTable();
double computeLik();
void evalBranchConfigs(vector<int> popsVec);
void readConfigFile();
void parseCmdLine(char* argv[]);
void checkConfigOptions();
double optimize_wrapper_nlopt(const vector<double> &vars, vector<double> &grad, void *data);
double check_constraints(const vector<double> &vars, vector<double> &grad, void *data);
void exploreProfiles (size_t &varIdx);


extern "C" {
int main_ms_ABLE(int ms_argc, char *ms_argv[], double ***onetreePoisTable);
double *** d3matrix(int x, int y, int z);
void freed3matrix(double ***m, int x, int y);
double ** d2matrix(int x, int y);
void freed2matrix(double **m, int x);
int ** d2int_matrix(int x, int y);
void freed2int_matrix(int **m, int x);

double ran1();
int getPopSampleStatus(int pop);
int getBrConfigNum(int *brConfVec);
}

#endif /* MAIN_H_ */

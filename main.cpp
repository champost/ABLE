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
map<vector<int>, int> intVec2BrConfig;
vector<int> npopVec;
map<vector<int>, double> dataConfigs;
map<vector<int>, double> selectConfigsMap;
map<long int, double> allConfigsMap;

int ms_argc;
char **ms_argv;
int ntrees;
int estimate;
long int finalTableSize;


double ranMT() { return(rMT()); }

//	conversion from decimal to base-mutClass
vector<int> getMutConfigVec(long int i) {
	int quo = i;
	int rem = 0;
	vector<int> MutConfig;
	for (int j = 0; j < brClass; j++) {
		if (quo) {
			rem = quo % mutClass;
			quo /= mutClass;
			MutConfig.push_back(rem);
		}
		else
			MutConfig.push_back(0);
	}
	reverse(MutConfig.begin(),MutConfig.end());

	return MutConfig;
}


//	conversion from decimal to base-mutClass
string getMutConfigStr(long int i) {
	int quo = i;
	int rem = 0;
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
		}
		else
			stst << "0";

		if (j < brClass - 1)
			stst << ",";
	}
	stst << "(";
	string config = stst.str();
	reverse(config.begin(),config.end());

	return config;
}


string getMutConfigStr(vector<int> configVec) {
	int size = configVec.size();
	stringstream stst;
	stst << "(";
	for (int j = 0; j < size-1; j++)
		stst << configVec[j] << ",";
	stst << configVec[size-1] << ")";

	return stst.str();
}


int getBrConfigNum(int *brConfVec) {
	vector<int> vec(brConfVec, brConfVec+npopVec.size());
	return intVec2BrConfig[vec];
}


double computeLik() {
	// calling ms
	main_ms(ms_argc, ms_argv);

	double loglik = 0.0;
	if (estimate == 2) {
		ofstream ifs("tmp.txt",ios::out);
		for (map<vector<int>, double>::iterator it = selectConfigsMap.begin(); it != selectConfigsMap.end(); it++) {
			ifs << getMutConfigStr(it->first) << " : " << it->second / ntrees << endl;
			loglik += log(it->second / ntrees) * dataConfigs[it->first];
		}
		ifs.close();
	}
	else if (estimate == 1) {
		for (map<vector<int>, double>::iterator it = selectConfigsMap.begin(); it != selectConfigsMap.end(); it++)
			loglik += log(it->second / ntrees) * dataConfigs[it->first];
	}
	else if (estimate == 0) {
		for (map<long int, double>::iterator it = allConfigsMap.begin(); it != allConfigsMap.end(); it++)
			printf("%s : %.5e\n", getMutConfigStr(it->first).c_str(), it->second/ntrees);
	}

	return loglik;
}


void calcFinalTable(double **onetreeTable) {
	if (estimate) {
		for (map<vector<int>, double>::iterator it = dataConfigs.begin(); it != dataConfigs.end(); it++) {
		    double jointPoisson = 1.0;

		    vector<int> vec = it->first;
			for (int j = 0; j < brClass; j++)
				jointPoisson *= onetreeTable[j][vec[j]];

			if (jointPoisson > 0.0)
				selectConfigsMap[vec] += jointPoisson;
		}

	}
	else {
		for (long int i = 0; i < finalTableSize; i++) {
		    double jointPoisson = 1.0;

		    vector<int> vec = getMutConfigVec(i);
			for (int j = 0; j < brClass; j++)
				jointPoisson *= onetreeTable[j][vec[j]];

			if (jointPoisson > 0.0)
				allConfigsMap[i] += jointPoisson;
		}
	}
}


double optimize_wrapper(const gsl_vector *vars, void *obj) {

	main_theta = gsl_vector_get(vars, 0);
	main_rho = gsl_vector_get(vars, 1);
	main_div_time = gsl_vector_get(vars, 2);
	if ((main_theta < 0) or (main_rho < 0) or (main_div_time < 0))
		return 999999;
	else
		return -computeLik();
}


void Tokenize(const string& str, vector<string>& tokens, const string& delimiters){
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}


void TrimSpaces(string& str)  {
	// Trim Both leading and trailing spaces
	size_t startpos = str.find_first_not_of(" \t\r\n"); // Find the first character position after excluding leading blank spaces
	size_t endpos = str.find_last_not_of(" \t\r\n"); // Find the first character position from reverse af

	// if all spaces or empty return an empty string
	if(( string::npos == startpos ) || ( string::npos == endpos))
	{
		str = "";
	}
	else
		str = str.substr( startpos, endpos-startpos+1 );
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


void readDataConfigs() {
	string line, del, keyVec;
	vector<string> tokens;
	vector<int> config;
	double val;

	ifstream ifs("dataconfig.txt",ios::in);
	while (getline(ifs,line)) {
		del = ":";
		tokens.clear();
		Tokenize(line, tokens, del);
		for(unsigned int j=0;j<tokens.size();j++){
			TrimSpaces(tokens[j]);
		}
		keyVec = tokens[0];
		val = atof(tokens[1].c_str());

		del = "(,)";
		tokens.clear();
		Tokenize(keyVec, tokens, del);
		for(unsigned int j=0;j<tokens.size();j++){
			TrimSpaces(tokens[j]);
			config.push_back(atoi(tokens[j].c_str()));
		}
		dataConfigs[config] = val;
		config.clear();
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
//		stringstream stst;
//		stst << ")";
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
//				stst << rem;
				sumConfig += rem;
				vec.push_back(rem);
			}
			else {
//				stst << "0";
				vec.push_back(0);
			}

//			if (j < npopVec.size() - 1)
//				stst << ",";
		}

		if (sumConfig == totPopSum)
			break;

		if (!skipConfig) {
//			stst << "(";
//			string config = stst.str();
//			reverse(config.begin(),config.end());
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
//	int nsam = atoi(argv[1]);
	int kmax = atoi(argv[argc-3]), npopSize = atoi(argv[argc-2]);
	char brFold = argv[argc-1][0];
	ms_argc = argc - 4;
	ms_argv = argv;
	ntrees = atoi(argv[2]);
	estimate = atoi(argv[argc-4]);

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

	if (estimate == 2) {
		readDataConfigs();

		main_theta = 1.0;
		main_rho = 1.0;
		main_div_time = 1.25;
		computeLik();


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
	else if (estimate == 1) {

		readDataConfigs();

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

		gsl_vector_set (x,0,ranMT());
		gsl_vector_set (x,1,ranMT());
		gsl_vector_set (x,2,ranMT());
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
			status = gsl_multimin_test_size (size, 0.001); //since we want more precision

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
//				printf ("LnL = %.6f size = %.6f\n", -s->fval, size);
////				printf ("LnL = %.6f\n", -s->fval);
//			}

		} while (status == GSL_CONTINUE && iter < 1000);

		gsl_vector_free(x);
		gsl_vector_free(ss);
		gsl_multimin_fminimizer_free (s);
	}
	else if (estimate == 0) {
		finalTableSize = (long int) pow(mutClass, brClass);
//		printf("\n\nfinalTableSize : %d", finalTableSize);

		time(&likStartTime);

		computeLik();

		time(&likEndTime);
		printf("\n\nTime taken for computation only : %.5f s", float(likEndTime - likStartTime));
	}

	return 0;
}




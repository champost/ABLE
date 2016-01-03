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

#include <nlopt.hpp>

#include "MersenneTwister.h"
#include "main.h"

using namespace std;

//************ EXTERN **************
int brClass, mutClass, foldBrClass = 0, allBrClasses;
//**********************************

MTRand rMT;
map<vector<int>, int> intVec2BrConfig;
vector<int> npopVec;
map<vector<int>, double> dataConfigs;
map<vector<int>, double> selectConfigsMap;
map<unsigned long int, double> allConfigsMap;
//vector<double> allConfigsLnL;
map<string, vector<int> > tbiIdx;
map<string, double> tbiStartVal;
map<string, int> parOrder;
map<double, vector<double> > bestGlobalSearchPointsMap, bestLocalSearchResultsMap;
string dataConfigFile;
int npops = 0, kmax = 0;
ofstream testLik, testConfig;

int ms_argc = 0;
char **ms_argv;
int ntrees = 0, treesSampled = 0;
int globalTrees = 500, localTrees = 1500, globalReps = 400, bestGlobalSearchPoints = 1;
double globalUpper = 5, globalLower = 1e-3;
bool skipGlobal = false, globalSearch = true, bSFS;
int estimate = 0, evalCount = 0;
unsigned long int finalTableSize;

double ranMT() { return(rMT()); }

//	conversion from decimal to base-mutClass
vector<int> getMutConfigVec(unsigned long int i) {

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
string getMutConfigStr(unsigned long int i) {

	int quo = i;
	int rem = 0;
	stringstream stst;
	stst << ")";
	for (int j = 0; j < brClass; j++) {
		if (quo) {
			rem = quo % mutClass;
			quo /= mutClass;
//			if (rem == mutClass-1)
//				stst << rem-1 << ">";
//			else
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

	treesSampled = 0;

	// calling ms
//	while ((double(selectConfigsMap.size()) < (0.8 * dataConfigs.size())) && (treesSampled < 10000))
		main_ms(ms_argc, ms_argv);

	double loglik = 0.0;
	if (bSFS) {
		ofstream ofs("bSFS.txt",ios::out);
		for (map<vector<int>, double>::iterator it = selectConfigsMap.begin(); it != selectConfigsMap.end(); it++) {
			ofs << getMutConfigStr(it->first) << " : " << scientific << it->second / treesSampled << endl;
			loglik += log(it->second / treesSampled) * dataConfigs[it->first];
		}
		ofs.close();


		loglik = loglik * dataConfigs.size() / selectConfigsMap.size();


		selectConfigsMap.clear();
	}
	else if (estimate > 0) {
		for (map<vector<int>, double>::iterator it = selectConfigsMap.begin(); it != selectConfigsMap.end(); it++)
			loglik += log(it->second / treesSampled) * dataConfigs[it->first];


		loglik = loglik * dataConfigs.size() / selectConfigsMap.size();


		selectConfigsMap.clear();
	}
	else if (estimate == 0) {
		ofstream ofs("expected_bSFS.txt",ios::out);
		for (map<unsigned long int, double>::iterator it = allConfigsMap.begin(); it != allConfigsMap.end(); it++)
			ofs << getMutConfigStr(it->first) << " : " << scientific << it->second / treesSampled << endl;
		ofs.close();

		allConfigsMap.clear();
	}

//	double average = 0.0, variance = 0.0;
//	for (int j = 0; j < treesSampled; j++)
//		average += allConfigsLnL[j];
//	average /= treesSampled;
//
//	for (int j = 0; j < treesSampled; j++)
//		variance += pow((average - allConfigsLnL[j]),2);
//	variance /= (treesSampled - 1);
//	printf("avgLnL : %.5f; varLnL : %.5f\n", average, variance);
//	allConfigsLnL.clear();

	return loglik;
}


void calcFinalTable(double **onetreeTable) {

	++treesSampled;

	double loglik = 0.0;

	if (bSFS || (estimate == 2)) {
		for (map<vector<int>, double>::iterator it = dataConfigs.begin(); it != dataConfigs.end(); it++) {
		    double jointPoisson = 1.0;

		    vector<int> vec = it->first;
			for (int j = 0; j < brClass; j++)
				jointPoisson *= onetreeTable[j][vec[j]];

			if (jointPoisson > 0.0)
				selectConfigsMap[vec] += jointPoisson;
		}
	}
	else if (estimate == 1) {
		for (map<vector<int>, double>::iterator it = dataConfigs.begin(); it != dataConfigs.end(); it++) {
		    double jointPoisson = 1.0;

		    vector<int> vec = it->first;
			for (int j = 0; j < brClass; j++)
				jointPoisson *= onetreeTable[j][vec[j]];

			if (jointPoisson > 0.0)
				selectConfigsMap[vec] += jointPoisson;

			if (selectConfigsMap.count(vec))
				loglik += log(selectConfigsMap[vec] / treesSampled) * it->second;
		}
		testLik << scientific << loglik << endl;
		testConfig << scientific << (double) selectConfigsMap.size()/dataConfigs.size() << endl;
	}
	else {
		for (unsigned long int i = 0; i < finalTableSize; i++) {
		    double jointPoisson = 1.0;

		    vector<int> vec = getMutConfigVec(i);
			for (int j = 0; j < brClass; j++)
				jointPoisson *= onetreeTable[j][vec[j]];

			if (jointPoisson > 0.0)
				allConfigsMap[i] += jointPoisson;
		}
	}
}


double optimize_wrapper_nlopt(const vector<double> &vars, vector<double> &grad, void *data) {

	++evalCount;
	if (!grad.empty()) {
		cerr << "Cannot proceed with ABLE" << endl;
		cerr << "Gradient based optimization not yet implemented..." << endl;
		exit(-1);
	}

	for (map<string, vector<int> >::iterator it = tbiIdx.begin(); it != tbiIdx.end(); it++) {
		stringstream stst;
		stst << vars[parOrder[it->first]];
		for (size_t i = 0; i < it->second.size(); i++)
			stst >> ms_argv[it->second[i]];
	}

	double loglik = computeLik();
	
	if (globalSearch)
		bestGlobalSearchPointsMap[-loglik] = vars;

	printf("%5d ", evalCount);
	for (size_t i = 0; i < vars.size(); i++)
		printf("%.5e ", vars[i]);
	printf(" Trees: %d ", treesSampled);
	printf(" LnL: %.6f\n", loglik);

	return loglik;

//	return computeLik();
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


void readConfigFile(int argc, char* argv[]) {

	string line, del, keyWord;
	vector<string> tokens;
	ifstream ifs("config.txt",ios::in);
	while (getline(ifs,line)) {
		del = " ";
		tokens.clear();
		Tokenize(line, tokens, del);
		for(unsigned int j=0;j<tokens.size();j++)
			TrimSpaces(tokens[j]);

		if ((tokens[0][0] != '#') && (line.size() > 0)) {
			if (tokens[0] == "pops") {
				stringstream stst1(tokens[1]);
				stst1 >> npops;
				for (int i = 0; i < npops; i++) {
					stringstream stst2(tokens[i+2]);
					int tmp;
					stst2 >> tmp;
					npopVec.push_back(tmp);
				}

			}
			else if (tokens[0] == "params") {
				double val;
				stringstream stst(tokens[2]);
				stst >> val;
				tbiStartVal[tokens[1]] = val;
			}
			else if (tokens[0] == "datafile") {
				dataConfigFile = tokens[1];
			}
			else if (tokens[0] == "estimate") {
				stringstream stst(tokens[1]);
				stst >> estimate;
			}
			else if (tokens[0] == "kmax") {
				stringstream stst(tokens[1]);
				stst >> kmax;
			}
			else if (tokens[0] == "folded")
				foldBrClass = 1;
			else if (tokens[0] == "global_trees") {
				stringstream stst(tokens[1]);
				stst >> globalTrees;
			}
			else if (tokens[0] == "local_trees") {
				stringstream stst(tokens[1]);
				stst >> localTrees;
			}
			else if (tokens[0] == "global_reps") {
				stringstream stst(tokens[1]);
				stst >> globalReps;
			}
			else if (tokens[0] == "global_upper") {
				stringstream stst(tokens[1]);
				stst >> globalUpper;
			}
			else if (tokens[0] == "global_lower") {
				stringstream stst(tokens[1]);
				stst >> globalLower;
			}
			else if (tokens[0] == "skip_global_search") {
				skipGlobal = true;
			}
			else if (tokens[0] == "bSFS") {
				bSFS = true;
			}
			else if (tokens[0] == "best_global_points") {
				stringstream stst(tokens[1]);
				stst >> bestGlobalSearchPoints;
			}
			else {
				cerr << "Unrecognised keyword \"" << tokens[0] << "\" found in the config file!" << endl;
				cerr << "Aborting ABLE..." << endl;
				exit(-1);
			}
		}
	}
	ifs.close();

	ms_argv = (char **)malloc( argc*sizeof(char *) ) ;
	for(int i =0; i < argc; i++)
		ms_argv[i] = (char *)malloc(30*sizeof(char) ) ;

	for (int i = 0; i < argc; i++) {
		string param(argv[i]);
		stringstream stst;
		if (param.substr(0,3) == "tbi") {
			tbiIdx[param].push_back(i);
			if (tbiStartVal.find(param) != tbiStartVal.end()) {
				stst << tbiStartVal[param];
				stst >> ms_argv[i];
			}
			else {
				double tmpPar = ranMT();
				tbiStartVal[param] = tmpPar;
				stst << tmpPar;
				stst >> ms_argv[i];
			}
		}
		else if ((i == 2) && (estimate > 0) && (globalTrees != 400)) {
			stst << globalTrees;
			stst >> ms_argv[2];
		}
		else {
			stst << argv[i];
			stst >> ms_argv[i];
		}
	}

//		for(int i = 1; i < ms_argc; i++)
//			printf("%s ",ms_argv[i]);
//		printf("\n");
//		exit(-1);

	if ((estimate == 2) && !tbiIdx.size()) {
		cerr << "Cannot proceed with inference" << endl;
		cerr << "\"tbi\" (To be Inferred) keywords need to be specified when \"estimate = 2\"" << endl;
		cerr << "Exiting ABLE..." << endl;
		exit(-1);
	}
}


void readDataConfigs() {
	string line, del, keyVec;
	vector<string> tokens;
	vector<int> config;
	double val;

	ifstream ifs(dataConfigFile.c_str(),ios::in);
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

	mutClass = kmax+2;
	brClass = npopVec[0]+1;
	for (size_t i = 1; i < npopVec.size(); i++)
		brClass *= (npopVec[i]+1);
	brClass -= 2;
	allBrClasses = brClass;

	//	fold the branch classes
	if (foldBrClass) {
		if (brClass % 2)
			brClass = (brClass+1) / 2;
		else
			brClass = brClass / 2;
	}

	for (size_t i = 1; i < npopVec.size(); i++) {
		totPopSum += npopVec[i];
		if (npopVec[i] > maxPopSize)
			maxPopSize = npopVec[i];
	}
	++maxPopSize;

	for (unsigned long int i = 1; i <= (unsigned long int) pow(maxPopSize,npopVec.size()); i++) {
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
	ntrees = atoi(argv[2]);
	ms_argc = argc;

	readConfigFile(argc, argv);

	evalBranchConfigs();

	if (estimate == 2) {

		readDataConfigs();

		double maxLnL;
		vector<double> parVec;
		int parCount = 0;
		for (map<string, double>::iterator it = tbiStartVal.begin(); it != tbiStartVal.end(); it++) {
			parVec.push_back(it->second);
			parOrder[it->first] = parCount;
			++parCount;
		}

		nlopt::opt opt(nlopt::GN_DIRECT_NOSCAL, tbiIdx.size());
		opt.set_lower_bounds(globalLower);
		opt.set_upper_bounds(globalUpper);
		opt.set_max_objective(optimize_wrapper_nlopt, NULL);
		opt.set_maxeval(globalReps);
		if (pow(4,parVec.size()) > globalReps)
			opt.set_maxeval(pow(4,parVec.size()));


		nlopt::opt local_opt(nlopt::LN_SBPLX, tbiIdx.size());
		local_opt.set_max_objective(optimize_wrapper_nlopt, NULL);
		local_opt.set_lower_bounds(globalLower);
		local_opt.set_upper_bounds(globalUpper);

		local_opt.set_maxeval(10);

//		local_opt.set_xtol_rel(1e-4);
		local_opt.set_initial_step(1);


		if (!skipGlobal) {
			printf("Starting global search: \n");

			evalCount = 0;
			time(&likStartTime);
			opt.optimize(parVec, maxLnL);

			printf("\nUsing the global search result(s) after %d evaluations as the starting point(s) for a refined local search...\n\n", evalCount);
			stringstream stst;
			stst << localTrees;
			stst >> ms_argv[2];

			int totEvalCount = 0;
			map<double, vector<double> >::iterator it = bestGlobalSearchPointsMap.begin();
			for (int localStart = 1; localStart <= bestGlobalSearchPoints; localStart++) {
				parVec = it->second;
				globalSearch = false;
				if (bestGlobalSearchPoints > 1)
					printf("Starting local search %d\n", localStart);

				evalCount = 0;
				local_opt.optimize(parVec, maxLnL);

				printf("Found the local maximum after %d evaluations\n", evalCount);
				printf("Found a maximum at ");
				for (size_t i = 0; i < parVec.size(); i++)
					printf("%.6f ", parVec[i]);
				printf("LnL = %.6f\n\n", maxLnL);

				bestLocalSearchResultsMap[-maxLnL] = parVec;
				totEvalCount += evalCount;
				++it;
			}
			if (bestGlobalSearchPoints > 1) {
				printf("\nFound the BEST local maximum after a total of %d evaluations\n", totEvalCount);
				printf("Found a maximum at ");
				it = bestLocalSearchResultsMap.begin();
				for (size_t i = 0; i < it->second.size(); i++)
					printf("%.6f ", it->second[i]);
				printf("LnL: %.6f \n\n", -it->first);
			}
		}
		else {
			printf("Skipping global search!\n");
			printf("\nUsing the user-specified/default values as a starting point for a local search...\n\n");
			time(&likStartTime);

			stringstream stst;
			stst << localTrees;
			stst >> ms_argv[2];

			globalSearch = false;
			evalCount = 0;
			local_opt.optimize(parVec, maxLnL);

			printf("Found the local maximum after %d evaluations\n", evalCount);
			printf("Found a maximum at ");
			for (size_t i = 0; i < parVec.size(); i++)
				printf("%.6f ", parVec[i]);
			printf("LnL = %.6f\n\n", maxLnL);
		}

		time(&likEndTime);
		printf("\nOverall time taken for optimization : %.5f s\n\n", float(likEndTime - likStartTime));
	}
	else if ((estimate == 1) || bSFS) {

		readDataConfigs();

		if (!tbiStartVal.empty()) {
			printf("Evaluating point likelihood at : \n");
			for (map<string, double>::iterator it = tbiStartVal.begin(); it != tbiStartVal.end(); it++)
				printf("%.6f ", it->second);
			printf("\n");
		}
		else {
			for (int i = 0; i < argc; i++)
				cout << argv[i] << " ";
			cout << endl;
		}

		if (!bSFS) {
			testLik.open("logliks.txt",ios::out);
			testConfig.open("propConfigs.txt",ios::out);
			time(&likStartTime);
			printf("LnL : %.6f\n", computeLik());
			time(&likEndTime);
			testLik.close();
			testConfig.close();
		}
		else {
			time(&likStartTime);
			printf("LnL : %.6f\n", computeLik());
			time(&likEndTime);
		}

		printf("Time taken for computation : %.5f s\n\n", float(likEndTime - likStartTime));
	}
	else if (estimate == 0) {

		finalTableSize = (long int) pow(mutClass, brClass);
		if (finalTableSize > 1000000000) { //	1 billion
			cerr << "\nThis is going to be too long a table to compute!" << endl;
			cerr << "Please contact the author Champak B. Reddy (champak.br@gmail.com) if you are really keen on going ahead with this" << endl;
			cerr << "Exiting ABLE...\n" << endl;
			exit(-1);
		}
//		printf("\n\nfinalTableSize : %.0f", finalTableSize);

		time(&likStartTime);
		computeLik();
		time(&likEndTime);

		printf("\n\nTime taken for computation : %.5f s\n", float(likEndTime - likStartTime));
	}

	return 0;
}




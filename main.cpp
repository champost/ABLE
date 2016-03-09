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
#include <numeric>
#include <limits>

#include <nlopt.hpp>
#include <gsl/gsl_rng.h>
#include <omp.h>

#include "main.h"
#include "utils.h"

using namespace std;

//************ EXTERN **************
int brClass, mutClass, foldBrClass = 0, allBrClasses, sampledPopsSize, ms_crash_flag;
//**********************************

gsl_rng * PRNG;

map<vector<int>, int> intVec2BrConfig;
map<int, vector<int> > tbiMsCmdIdx;
map<int, double> tbiUserVal;
map<int, vector<double> > tbiSearchBounds;
map<int, int> trackSelectConfigs, parConstraints;

string dataConfigFile, configFile;
ofstream testLik, testConfig;

vector<int> allConfigs, trackSelectConfigsForInf, sampledPops, allPops;
vector<vector<int> > dataConfigs;
vector<double> dataConfigFreqs, selectConfigFreqs, allConfigFreqs, upperBounds, lowerBounds, bestGlobalSPars, lastValidPars;
vector<gsl_rng *> PRNGThreadVec;

nlopt::opt opt;
nlopt::opt local_opt;

char **ms_argv;

double **onetreePoisTable;
#pragma omp threadprivate(onetreePoisTable)

int ms_argc = 0;
int npops = 0, kmax = 0;
int estimate = 0, evalCount = 0;
int globalTrees = 0, localTrees = 0, globalEvals = 0, localEvals = 0, refineLikTrees = 0, ms_trees = 1, reportEveryEvals = 0;

double globalUpper = 5, globalLower = 1e-3, penLnL, dataLnL, bestGlobalSlLnL, lastValidLnL;
bool skipGlobal = false, bSFS = false, profileLikBool = true, onlyProfiles = false, abortNLopt = false, seedPRNGBool = false;
unsigned long int finalTableSize, seedPRNG;

enum SearchStates {GLOBAL, LOCAL, OTHER};
SearchStates currState = OTHER;


double ran1()
{
	return gsl_rng_uniform(PRNGThreadVec[omp_get_thread_num()]);
}

void free_objects() {
	for(int i = 0; i < ms_argc; i++)
		free(ms_argv[i]);
	free(ms_argv);

	gsl_rng_free(PRNG);
	for (size_t i = 0; i < PRNGThreadVec.size(); i++)
		gsl_rng_free(PRNGThreadVec[i]);

}

void profileLik(vector<double> MLEparVec) {

	double parLowerBound, parUpperBound, MLEparVal, MLEparLik;

	stringstream ststTree;
	ststTree << localTrees;
	ststTree >> ms_argv[2];

	size_t parIdx = 0;
	//	Likelihood profiles for all parameters to be inferred (tbi)
	for (map<int, vector<int> >::iterator it = tbiMsCmdIdx.begin(); it != tbiMsCmdIdx.end(); it++) {
		stringstream tbiFileName;
		double quarterGlobalRange = (upperBounds[parIdx] - lowerBounds[parIdx]) / 4;
		MLEparVal = MLEparVec[it->first-1];
		parUpperBound = MLEparVal + quarterGlobalRange;
		parLowerBound = MLEparVal - quarterGlobalRange;
		if (parLowerBound < 0)
			parLowerBound = MLEparVal / 2;

		tbiFileName << "tbi" << it->first;
		printf("\nCalculating profiles of the likelihood surface for %s\n", tbiFileName.str().c_str());

		tbiFileName << ".txt";
		ofstream outFile(tbiFileName.str().c_str(),ios::out);

		//	Calculating the LnL at the MLE
		for (size_t i = 0; i < it->second.size(); i++) {
			stringstream ststMLEPar;
			ststMLEPar << MLEparVal;
			ststMLEPar >> ms_argv[it->second[i]];
		}
		MLEparLik = computeLik();
		outFile << scientific << MLEparVal << "\t" << MLEparLik << endl;
		bool insertMLEpar = true;

		//	Likelihood profiles for this parameter
		vector<double> parRange = logspaced(parLowerBound, parUpperBound, 10);
		for (size_t j = 0; j < parRange.size(); j++) {
			for (size_t i = 0; i < it->second.size(); i++) {
				stringstream ststProfilePar;
				ststProfilePar << parRange[j];
				ststProfilePar >> ms_argv[it->second[i]];
			}
			if (insertMLEpar && (parRange[j] > MLEparVal)) {
				outFile << scientific << MLEparVal << "\t" << MLEparLik << endl;
				insertMLEpar = false;
			}
			outFile << scientific << parRange[j] << "\t" << computeLik() << endl;
		}

		outFile.close();

		//	Resetting the MLE for this parameter
		for (size_t i = 0; i < it->second.size(); i++) {
			stringstream ststMLEPar;
			ststMLEPar << MLEparVal;
			ststMLEPar >> ms_argv[it->second[i]];
		}
		++parIdx;
	}
}


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
/*
			if (rem == mutClass-1)
				stst << rem-1 << ">";
			else
*/
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


int getPopSampleStatus(int pop) { return allPops[pop]; }


int getBrConfigNum(int *brConfVec) {

	if (accumulate(brConfVec, brConfVec+sampledPops.size(), 0)) {
		vector<int> vec(brConfVec, brConfVec+sampledPops.size());
		return intVec2BrConfig[vec];
	}

	return -1;
}


void calcBSFSTable() {

	ms_crash_flag = 0;
	{
		stringstream stst;
		stst << 1;
		stst >> ms_argv[2];
	}

	if (bSFS || (estimate == 2)) {
#pragma omp parallel for
		for (int trees = 0; trees < ms_trees; trees++) {
			if (!ms_crash_flag) {
				onetreePoisTable = d2matrix(brClass, mutClass);
				// calling ms for sampling genealogies
				main_ms_ABLE(ms_argc, ms_argv, onetreePoisTable);

				for (size_t i = 0; i < dataConfigs.size(); i++) {
				    double jointPoisson = 1.0;

					for (int j = 0; j < brClass; j++)
						jointPoisson *= onetreePoisTable[j][dataConfigs[i][j]];

					if (jointPoisson > (numeric_limits<double>::min()*ms_trees)) {
						selectConfigFreqs[i] += jointPoisson;
						trackSelectConfigsForInf[i] = 1;
					}
				}
				freed2matrix(onetreePoisTable, brClass);
			}
		}
	}
	else if (estimate == 1) {
		int treesSampled = 0;
#pragma omp parallel for
		for (int trees = 0; trees < ms_trees; trees++) {
			if (!ms_crash_flag) {
				onetreePoisTable = d2matrix(brClass, mutClass);
				// calling ms for sampling genealogies
				main_ms_ABLE(ms_argc, ms_argv, onetreePoisTable);

#pragma omp critical
				{
					++treesSampled;
					double loglik = 0.0;
					for (size_t i = 0; i < dataConfigs.size(); i++) {
						double jointPoisson = 1.0;

						for (int j = 0; j < brClass; j++)
							jointPoisson *= onetreePoisTable[j][dataConfigs[i][j]];

						if (jointPoisson > (numeric_limits<double>::min()*ms_trees))
							selectConfigFreqs[i] += jointPoisson;

						if (selectConfigFreqs[i] != 0.0) {
							loglik += log(selectConfigFreqs[i] / treesSampled) * dataConfigFreqs[i];
							trackSelectConfigs[i] = 1;
						}
					}
					testLik << scientific << loglik << endl;
					testConfig << scientific << (double) trackSelectConfigs.size()/dataConfigs.size() << endl;
				}
				freed2matrix(onetreePoisTable, brClass);
			}
		}
	}
	else {
		allConfigs = vector<int>(finalTableSize,-1);
		allConfigFreqs = vector<double>(finalTableSize,0.0);

#pragma omp parallel for
		for (int trees = 0; trees < ms_trees; trees++) {
			if (!ms_crash_flag) {
				onetreePoisTable = d2matrix(brClass, mutClass);
				// calling ms for sampling genealogies
				main_ms_ABLE(ms_argc, ms_argv, onetreePoisTable);

				for (unsigned long int i = 0; i < finalTableSize; i++) {
				    double jointPoisson = 1.0;

				    vector<int> vec = getMutConfigVec(i);
					for (int j = 0; j < brClass; j++)
						jointPoisson *= onetreePoisTable[j][vec[j]];

					if (jointPoisson > (numeric_limits<double>::min()*ms_trees)) {
						allConfigs[i] = i;
						allConfigFreqs[i] += jointPoisson;
					}
				}
				freed2matrix(onetreePoisTable, brClass);
			}
		}
	}
}


double computeLik() {

	double loglik = 0.0;

	selectConfigFreqs = vector<double>(dataConfigFreqs.size(),0.0);
	trackSelectConfigsForInf = vector<int>(dataConfigFreqs.size(),0);

	//	calculating the bSFS config. probs.
	calcBSFSTable();

	if (!ms_crash_flag) {

		int trackedConfigs = 0;
		if (bSFS || (estimate > 1))
			trackedConfigs = accumulate(trackSelectConfigsForInf.begin(),trackSelectConfigsForInf.end(),0);
		else if (estimate == 1)
			trackedConfigs = trackSelectConfigs.size();

		if (bSFS) {
			ofstream ofs("bSFS.txt",ios::out);
			for (size_t i = 0; i < dataConfigs.size(); i++) {
				if (selectConfigFreqs[i] != 0.0) {
					ofs << getMutConfigStr(dataConfigs[i]) << " : " << scientific << selectConfigFreqs[i] / ms_trees << endl;
					loglik += log(selectConfigFreqs[i] / ms_trees) * dataConfigFreqs[i];
				}
			}
			ofs.close();

			if (trackedConfigs)
				loglik *= (double) dataConfigFreqs.size() / trackedConfigs;
		}
		else if (estimate > 0) {
			for (size_t i = 0; i < dataConfigs.size(); i++) {
				if (selectConfigFreqs[i] != 0.0)
					loglik += log(selectConfigFreqs[i] / ms_trees) * dataConfigFreqs[i];
			}

			if (trackedConfigs)
				loglik *= (double) dataConfigFreqs.size() / trackedConfigs;
		}
		else if (estimate == 0) {
			ofstream ofs("expected_bSFS.txt",ios::out);
			for (size_t i = 0; i < allConfigs.size(); i++)
				if (allConfigs[i] >= 0)
					ofs << getMutConfigStr(allConfigs[i]) << " : " << scientific << allConfigFreqs[allConfigs[i]] / ms_trees << endl;
			ofs.close();

			allConfigs.clear();
			allConfigFreqs.clear();
		}
	}

	trackSelectConfigsForInf.clear();
	trackSelectConfigs.clear();
	selectConfigFreqs.clear();

	return loglik;
}


double optimize_wrapper_nlopt(const vector<double> &vars, vector<double> &grad, void *data) {

	bool parConstraintPass = true;
	double loglik = 0.0;

	if (!grad.empty()) {
		cerr << "Cannot proceed with ABLE" << endl;
		cerr << "Gradient based optimization not yet implemented..." << endl;
		free_objects();
		exit(-1);
	}

	//	global search exit status : 1234567
	if ((currState == GLOBAL) && (evalCount > globalEvals))
			return 1234567;

	//	checking for simple non linear constraints between the free params
	for (map<int , int>::iterator it = parConstraints.begin(); it != parConstraints.end(); it++) {
		if (vars[it->first-1] >= vars[it->second-1]) {
			parConstraintPass = false;
			break;
		}
	}

	//	adaptive penalization based on the difference from last sampled and valid parameter vector
	//	in order to better explore non-linearly constrained space
	double penFactor = 0.0;
	if (lastValidPars.size() > 0) {
		for (size_t i = 0; i < vars.size(); i++)
			penFactor += fabs(vars[i] - lastValidPars[i]);
	}

	//	the standard global/local LnL search procedure (when all else is OK)
	if (parConstraintPass) {
		for (map<int, vector<int> >::iterator it = tbiMsCmdIdx.begin(); it != tbiMsCmdIdx.end(); it++) {
			for (size_t i = 0; i < it->second.size(); i++) {
				stringstream stst;
				stst << vars[it->first-1];
				stst >> ms_argv[it->second[i]];
			}
		}

/*
		for(int i = 1; i < ms_argc; i++)
			printf("%s ",ms_argv[i]);
		printf("\n");
//		free_objects(); exit(-1);
*/

		loglik = computeLik();

		if (ms_crash_flag) {
/*
			printf("%5d ", evalCount);
			for (size_t i = 0; i < vars.size(); i++)
				printf("%.5e ", vars[i]);
*/

			if (penLnL == 0)
				penLnL = 100*dataLnL;
			else {
				if (lastValidPars.size() > 0)
					penLnL = lastValidLnL + penFactor*dataLnL;
				else
					penLnL += dataLnL;
			}

/*
			printf(" Trees: %d ", 0);
			printf(" Penalised LnL: %.6f ", penLnL);
			printf(" ms CRASH!\n");
*/

			return penLnL;
		}

		++evalCount;
		printf("%5d ", evalCount);
		for (size_t i = 0; i < vars.size(); i++)
			printf("%.5e ", vars[i]);
		printf(" Trees: %d ", ms_trees);
		printf(" LnL: %.6f\n", loglik);

		if (loglik == 0.0) {
			//	global/local search exit status : 1234567
			abortNLopt = true;
			if (currState == GLOBAL)
				return 1234567;
			else if (currState == LOCAL)
				return 7654321;
		}

		if ((currState == GLOBAL) && reportEveryEvals && !(evalCount % reportEveryEvals) && (evalCount < globalEvals)) {
			printf("\nReporting the best MLE after %d evaluations\n", evalCount);
			for (size_t i = 0; i < bestGlobalSPars.size(); i++)
				printf("%.6f ", bestGlobalSPars[i]);
			printf("LnL = %.6f\n\n", bestGlobalSlLnL);
		}
	}
	//	penalising likelihood evaluation when searching outside the constrained zone
	//	special case of first search point is in the forbidden zone : 100*dataLnL
	else {
/*
		printf("%5d ", evalCount);
		for (size_t i = 0; i < vars.size(); i++)
			printf("%.5e ", vars[i]);
*/

		if (penLnL == 0)
			penLnL = 100*dataLnL;
		else {
			if (lastValidPars.size() > 0)
				penLnL = lastValidLnL + penFactor*dataLnL;
			else
				penLnL += dataLnL;
		}

/*
		printf(" Trees: %d ", 0);
		printf(" Penalised LnL: %.6f\n", penLnL);
*/

		return penLnL;
	}

	//	side stepping nlopt by storing the best LnL and parameters
	if ((currState == GLOBAL) && (loglik > bestGlobalSlLnL)) {
			bestGlobalSlLnL = loglik;
			bestGlobalSPars = vars;
	}

	lastValidLnL = loglik;
	lastValidPars = vars;

	return loglik;
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
		dataConfigs.push_back(config);
		dataConfigFreqs.push_back(val);

		config.clear();
	}
	ifs.close();

	dataLnL = 0.0;
	for (size_t i = 0; i < dataConfigs.size(); i++)
		dataLnL += log(dataConfigFreqs[i]) * dataConfigFreqs[i];

	bestGlobalSlLnL = 100000*dataLnL;
}


//	conversion from decimal to base-(maxPopSize+1)
void evalBranchConfigs() {

	int quo, rem, maxPopSize, totPopSum, count = 0, sumConfig;
	bool skipConfig;
	maxPopSize = totPopSum = sampledPops[0];

	mutClass = kmax+2;
	brClass = sampledPops[0]+1;
	for (size_t i = 1; i < sampledPops.size(); i++)
		brClass *= (sampledPops[i]+1);
	brClass -= 2;
	allBrClasses = brClass;

	//	fold the branch classes
	if (foldBrClass) {
		if (brClass % 2)
			brClass = (brClass+1) / 2;
		else
			brClass = brClass / 2;
	}

	for (size_t i = 1; i < sampledPops.size(); i++) {
		totPopSum += sampledPops[i];
		if (sampledPops[i] > maxPopSize)
			maxPopSize = sampledPops[i];
	}
	++maxPopSize;

	for (unsigned long int i = 1; i <= (unsigned long int) pow(maxPopSize,sampledPops.size()); i++) {
		quo = i;
		rem = 0;
/*
		stringstream stst;
		stst << ")";
*/
		sumConfig = 0;
		skipConfig = false;
		vector<int> vec;

		for (size_t j = 0; j < sampledPops.size(); j++) {
			if (quo) {
				rem = quo % (maxPopSize);
				quo /= (maxPopSize);
				if (rem > sampledPops[sampledPops.size()-1-j]) {
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

/*
			if (j < sampledPops.size() - 1)
				stst << ",";
*/
		}

		if (sumConfig == totPopSum)
			break;

		if (!skipConfig) {
/*
			stst << "(";
			string config = stst.str();
			reverse(config.begin(),config.end());
*/
			reverse(vec.begin(),vec.end());
			intVec2BrConfig[vec] = count;
			++count;
/*
			printf("%d\t%d\t%s\n", i, count, config.c_str());
			printf("%d\t%s\n", count, config.c_str());
*/
		}
	}
}


void readConfigFile(char* argv[]) {

	string line, del, keyWord;
	vector<string> tokens;
	ifstream ifs(configFile.c_str(),ios::in);
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
					if (stst2.str() != "u") {
						int tmp;
						stst2 >> tmp;
						sampledPops.push_back(tmp);
						allPops.push_back(1);
					}
					else
						allPops.push_back(0);
				}
				sampledPopsSize = sampledPops.size();
			}
			else if (tokens[0] == "start") {
				double val;
				if (tokens[1] == "all") {
					for(unsigned int j = 2; j < tokens.size(); j++) {
						stringstream stst_val(tokens[j]);
						stst_val >> val;
						tbiUserVal[j-1] = val;
					}
				}
				else {
					stringstream stst(tokens[2]);
					stst >> val;
					tbiUserVal[atoi(tokens[1].substr(3).c_str())] = val;
				}
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
			else if (tokens[0] == "global_search_trees") {
				stringstream stst(tokens[1]);
				stst >> globalTrees;
			}
			else if (tokens[0] == "local_search_trees") {
				stringstream stst(tokens[1]);
				stst >> localTrees;
			}
			else if (tokens[0] == "global_search_evals") {
				stringstream stst(tokens[1]);
				stst >> globalEvals;
			}
			else if (tokens[0] == "local_search_evals") {
				stringstream stst(tokens[1]);
				stst >> localEvals;
			}
			else if (tokens[0] == "global_upper_bound") {
				stringstream stst(tokens[1]);
				stst >> globalUpper;
			}
			else if (tokens[0] == "global_lower_bound") {
				stringstream stst(tokens[1]);
				stst >> globalLower;
			}
			else if (tokens[0] == "skip_global_search") {
				skipGlobal = true;
			}
			else if (tokens[0] == "bounds") {
				int paramID = atoi(tokens[1].substr(3).c_str());
				for(unsigned int j = 2; j < 4; j++) {
					double val;
					stringstream stst(tokens[j]);
					stst >> val;
					tbiSearchBounds[paramID].push_back(val);
				}
			}
			else if (tokens[0] == "constrain") {
				int paramID1 = atoi(tokens[1].substr(3).c_str()), paramID2 = atoi(tokens[2].substr(3).c_str());
				parConstraints[paramID1] = paramID2;
			}
			else if (tokens[0] == "bSFS") {
				bSFS = true;
			}
			else if (tokens[0] == "no_profiles") {
				profileLikBool = false;
			}
			else if (tokens[0] == "only_profiles") {
				onlyProfiles = true;
			}
			else if (tokens[0] == "seed_PRNG") {
				stringstream stst(tokens[1]);
				stst >> seedPRNG;
				seedPRNGBool = true;
			}
			else if (tokens[0] == "refine_likelihoods") {
				stringstream stst(tokens[1]);
				stst >> refineLikTrees;
			}
			else if (tokens[0] == "report_likelihoods") {
				stringstream stst(tokens[1]);
				stst >> reportEveryEvals;
			}
			else {
				cerr << "Unrecognised keyword \"" << tokens[0] << "\" found in the config file!" << endl;
				cerr << "Aborting ABLE..." << endl;
				exit(-1);
			}
		}
	}
	ifs.close();

	ms_argv = (char **)malloc( ms_argc*sizeof(char *) ) ;
	for(int i =0; i < ms_argc; i++)
		ms_argv[i] = (char *)malloc(30*sizeof(char) ) ;

	for (int i = 0; i < ms_argc; i++) {
		string param(argv[i]);
		stringstream stst;
		if (param.substr(0,3) == "tbi") {
			int paramID = atoi(param.substr(3).c_str());
			tbiMsCmdIdx[paramID].push_back(i);
			if (tbiUserVal.find(paramID) != tbiUserVal.end()) {
				stst << tbiUserVal[paramID];
				stst >> ms_argv[i];
			}
			else {
				if (estimate == 2 && (onlyProfiles || skipGlobal)) {
					if (onlyProfiles)
						cerr << "\nCannot proceed with plotting only profiles" << endl;
					else if (skipGlobal)
						cerr << "\nCannot proceed with local search" << endl;
					cerr << "You need to specify values for the \"tbi\" keywords using the \"start\" keyword" << endl;
					cerr << "Exiting ABLE...\n" << endl;
					free_objects();
					exit(-1);
				}

				double tmpPar = gsl_rng_uniform(PRNG);
				tbiUserVal[paramID] = tmpPar;
				stst << tmpPar;
				stst >> ms_argv[i];
			}
		}
		else {
			stst << argv[i];
			stst >> ms_argv[i];
		}
	}

/*
	for(int i = 1; i < ms_argc; i++)
		printf("%s ",ms_argv[i]);
	printf("\n");
	free_objects();
	exit(-1);
*/

	if (estimate == 2) {
		if (!tbiMsCmdIdx.size()) {
			cerr << "\nCannot proceed with inference" << endl;
			cerr << "\"tbi\" (To be Inferred) keywords need to be specified when \"estimate = 2\"" << endl;
			cerr << "Exiting ABLE...\n" << endl;
			free_objects();
			exit(-1);
		}
		else {
			stringstream stst;
			if (globalTrees == 0)
				globalTrees = 1000*tbiMsCmdIdx.size();

			if (localTrees == 0)
				localTrees = 1000*tbiMsCmdIdx.size();
		}
	}

	if (onlyProfiles)
		skipGlobal = true;
}


int main(int argc, char* argv[]) {

	omp_set_num_threads(omp_get_num_procs());

	seedPRNG = hash(time(NULL), clock());
	for (int i = 0; i < omp_get_num_procs(); i++) {
		PRNGThreadVec.push_back(gsl_rng_alloc(gsl_rng_mt19937));
		gsl_rng_set(PRNGThreadVec[i], seedPRNG + i);
	}

	PRNG = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(PRNG, hash(time(NULL), clock()));

	time_t likStartTime, likEndTime;
/*
	int nsam = atoi(argv[1]);
	ntrees = atoi(argv[2]);
*/
	ms_argc = argc;
	configFile = string(argv[argc - 1]);
	if (configFile == "-T")
		configFile = "config.txt";
	else
		ms_argc = argc - 1;

	readConfigFile(argv);

	if (seedPRNGBool) {
		for (size_t i = 0; i < PRNGThreadVec.size(); i++)
			gsl_rng_set(PRNGThreadVec[i], seedPRNG + i);
	}

	evalBranchConfigs();

	if (estimate == 2) {

		readDataConfigs();

		double maxLnL;
		vector<double> parVec;
		for (map<int, double>::iterator it = tbiUserVal.begin(); it != tbiUserVal.end(); it++) {
			parVec.push_back(it->second);

			if (tbiSearchBounds.find(it->first) != tbiSearchBounds.end()) {
				lowerBounds.push_back(tbiSearchBounds[it->first][0]);
				upperBounds.push_back(tbiSearchBounds[it->first][1]);
			}
			else {
				lowerBounds.push_back(globalLower);
				upperBounds.push_back(globalUpper);
			}
		}

		opt = nlopt::opt(nlopt::GN_DIRECT, tbiMsCmdIdx.size());
		local_opt = nlopt::opt(nlopt::LN_SBPLX, tbiMsCmdIdx.size());

		opt.set_stopval(1234567);
		opt.set_lower_bounds(lowerBounds);
		opt.set_upper_bounds(upperBounds);
		opt.set_max_objective(optimize_wrapper_nlopt, NULL);
		if (!skipGlobal && (globalEvals < 3000*(int)tbiMsCmdIdx.size())) {
			if (globalEvals)
				printf("\nToo few global_search_points for the specified number of free parameters\nPlease consider increasing \"global_search_evals\"...\n");
			else
				globalEvals = 3000*tbiMsCmdIdx.size();
		}

		local_opt.set_stopval(7654321);
		local_opt.set_max_objective(optimize_wrapper_nlopt, NULL);
		local_opt.set_lower_bounds(lowerBounds);
		local_opt.set_upper_bounds(upperBounds);
		if (localEvals)
			local_opt.set_maxeval(localEvals);
		else
			local_opt.set_maxeval(10000);
		local_opt.set_xtol_rel(1e-2);
		local_opt.set_initial_step((globalUpper-globalLower)/5);

		if (!skipGlobal) {
			printf("\nStarting global search: \n");

			ms_trees = globalTrees;
			evalCount = 0;
			penLnL = 0;
			currState = GLOBAL;
			time(&likStartTime);

			try {
				opt.optimize(parVec, maxLnL);
				throw abortNLopt;
			}
			catch (...) {
				if (abortNLopt) {
					cerr << "Something went wrong in the calculation of the LnL!" << endl;
					cerr << "Aborting ABLE..." << endl;

					for (size_t i = 0; i < parVec.size(); i++)
						printf("%.6f ", parVec[i]);
					printf("LnL = %.6f\n\n", maxLnL);

					free_objects();
					exit(-1);
				}
			}

			printf("\nUsing the global search result(s) after %d evaluations as the starting point for a refined local search...\n\n", evalCount);

			ms_trees = localTrees;
			evalCount = 0;
			penLnL = 0;
			currState = LOCAL;
			parVec = bestGlobalSPars;

			try {
				local_opt.optimize(parVec, maxLnL);
				throw abortNLopt;
			}
			catch (...) {
				if (abortNLopt) {
					cerr << "Something went wrong in the calculation of the LnL!" << endl;
					cerr << "Aborting ABLE..." << endl;

					for (size_t i = 0; i < parVec.size(); i++)
						printf("%.6f ", parVec[i]);
					printf("LnL = %.6f\n\n", maxLnL);

					free_objects();
					exit(-1);
				}
			}

			printf("Found the local maximum after %d evaluations\n", evalCount);
			if (refineLikTrees) {
				printf("Refining the likelihood at the MLE using %d genealogies...\n", refineLikTrees);
				for (map<int, vector<int> >::iterator it = tbiMsCmdIdx.begin(); it != tbiMsCmdIdx.end(); it++) {
					for (size_t i = 0; i < it->second.size(); i++) {
						stringstream stst;
						stst << parVec[it->first-1];
						stst >> ms_argv[it->second[i]];
					}
				}

				ms_trees = refineLikTrees;
				maxLnL = computeLik();
			}

			printf("Found a maximum at ");
			for (size_t i = 0; i < parVec.size(); i++)
				printf("%.6f ", parVec[i]);
			printf("LnL = %.6f\n\n", maxLnL);

			time(&likEndTime);
			printf("\nOverall time taken for optimization : %.5f s\n\n", float(likEndTime - likStartTime));
		}
		else if (!onlyProfiles) {
			printf("Skipping global search!\n");
			printf("\nUsing the user-specified/default values as a starting point for a local search...\n\n");
			time(&likStartTime);

			ms_trees = localTrees;
			evalCount = 0;
			penLnL = 0;
			currState = LOCAL;

			try {
				local_opt.optimize(parVec, maxLnL);
				throw abortNLopt;
			}
			catch (...) {
				if (abortNLopt) {
					cerr << "Something went wrong in the calculation of the LnL!" << endl;
					cerr << "Aborting ABLE..." << endl;

					for (size_t i = 0; i < parVec.size(); i++)
						printf("%.6f ", parVec[i]);
					printf("LnL = %.6f\n\n", maxLnL);

					free_objects();
					exit(-1);
				}
			}

			printf("Found the local maximum after %d evaluations\n", evalCount);
			if (refineLikTrees) {
				printf("Refining the likelihood at the MLE using %d genealogies...\n", refineLikTrees);
				for (map<int, vector<int> >::iterator it = tbiMsCmdIdx.begin(); it != tbiMsCmdIdx.end(); it++) {
					for (size_t i = 0; i < it->second.size(); i++) {
						stringstream stst;
						stst << parVec[it->first-1];
						stst >> ms_argv[it->second[i]];
					}
				}

				ms_trees = refineLikTrees;
				maxLnL = computeLik();
			}

			printf("Found a maximum at ");
			for (size_t i = 0; i < parVec.size(); i++)
				printf("%.6f ", parVec[i]);
			printf("LnL = %.6f\n\n", maxLnL);

			time(&likEndTime);
			printf("\nOverall time taken for optimization : %.5f s\n\n", float(likEndTime - likStartTime));
		}

		// temporarily deactivated!
//		currState = OTHER;
//		if (onlyProfiles || profileLikBool)
//			profileLik(parVec);
	}
	else if ((estimate == 1) || bSFS) {

		readDataConfigs();

		printf("Evaluating point likelihood at : \n");
		if (!tbiUserVal.empty()) {
			for (map<int, double>::iterator it = tbiUserVal.begin(); it != tbiUserVal.end(); it++)
				printf("%.6f ", it->second);
			printf("\n");
		}
		else {
			for (int i = 0; i < ms_argc; i++)
				cout << argv[i] << " ";
			cout << endl;
		}

		{
			stringstream stst;
			stst << ms_argv[2];
			stst >> ms_trees;
		}

		currState = OTHER;
		double loglik;
		if (!bSFS) {
			testLik.open("logliks.txt",ios::out);
			testConfig.open("propConfigs.txt",ios::out);
			time(&likStartTime);
			loglik = computeLik();
			time(&likEndTime);
			testLik.close();
			testConfig.close();
		}
		else {
			time(&likStartTime);
			loglik = computeLik();
			time(&likEndTime);
		}

		if (ms_crash_flag) {
			cerr << "\nABLE failed to simulate genealogies with the demographic parameters that have been specified!" << endl;
			cerr << "Please consider changing them or contact the author Champak B. Reddy (champak.br@gmail.com)" << endl;
			cerr << "Exiting ABLE...\n" << endl;

			free_objects();
			exit(-1);
		}
		else {
			printf("LnL : %.6f (Trees sampled : %d)\n", loglik, ms_trees);
			printf("Time taken for computation : %.5f s\n\n", float(likEndTime - likStartTime));
		}
	}
	else if (estimate == 0) {

		finalTableSize = (long int) pow(mutClass, brClass);
		if (finalTableSize > 1000000000) { //	1 billion
			cerr << "\nThis is going to be too long a table to compute!" << endl;
			cerr << "Please contact the author Champak B. Reddy (champak.br@gmail.com) if you are really keen on going ahead with this" << endl;
			cerr << "Exiting ABLE...\n" << endl;
			free_objects();
			exit(-1);
		}
//		printf("\n\nfinalTableSize : %.0f", finalTableSize);

		{
			stringstream stst;
			stst << ms_argv[2];
			stst >> ms_trees;
		}
		currState = OTHER;
		time(&likStartTime);
		computeLik();
		time(&likEndTime);

		if (ms_crash_flag) {
			cerr << "\nABLE failed to simulate genealogies with the demographic parameters that have been specified!" << endl;
			cerr << "Please consider changing them or contact the author Champak B. Reddy (champak.br@gmail.com)" << endl;
			cerr << "Exiting ABLE...\n" << endl;

			free_objects();
			exit(-1);
		}
		else {
			printf("Trees sampled : %d\n", ms_trees);
			printf("\nTime taken for computation : %.5f s\n", float(likEndTime - likStartTime));
		}
	}

	free_objects();

	return 0;
}




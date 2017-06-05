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
#include <iomanip>
#include <set>

#include <nlopt.hpp>
#include <gsl/gsl_rng.h>
#include <omp.h>

#include "main.h"
#include "utils.h"

using namespace std;

//************ EXTERN **************
int brClass = 0, mutClass = 0, foldBrClass = 0, allBrClasses = 0, sampledPopsSize = 0, poisTableSize = 0;
int *poisTableScaleOrder;
int **blockLengthsMat;
//**********************************

gsl_rng * PRNG;

map<vector<int>, int> intVec2BrConfig;
map<int, vector<int> > tbiMsCmdIdx;
map<int, double> tbiUserVal;
map<int, vector<double> > tbiSearchBounds, tbiProfilesGrid;
map<int, int> parConstraints, tbi2ParVec, ParVecToTbi;
map<int, bool> setRandomPars;
map<double, vector<double> > bestParsMap;
//map<int, int> trackSelectConfigs;



//------------------------------------------------------
vector<string> dataFile;
int mbSFSLen = 0;
vector<int> dataScales;
set<int> uniqueScales;
vector<vector<int> > data2PoisTableIdx;



vector<vector<vector<int> > > dataConfigs;
vector<vector<double> > dataConfigFreqs, selectConfigFreqs;
vector<vector<int> > trackSelectConfigsForInf;
//------------------------------------------------------



string dataFileFormat = "bSFS", alleleType = "genotype", configFile, globalSearchAlg, bSFSFile = "bSFS.txt", data2bSFSFile;
//ofstream testLik, testConfig;

vector<int> allConfigs, sampledPops, allPops, profileVarKey;
vector<double> allConfigFreqs, upperBounds, lowerBounds, hardUpperBounds, hardLowerBounds, bestGlobalSPars, bestLocalSPars, profileVars;
vector<string> cmdLine;
vector<gsl_rng *> PRNGThreadVec;

nlopt::opt opt, local_opt, AUGLAG;

char **ms_argv;

int ms_argc = 0, ms_crash_flag = 0;
int npops = 0, kmax = 0;
int estimate = 0, evalCount = 0, crash_counter = 0, sampledTrees = 0, recLen = 0;
int globalTrees = 0, localTrees = 0, globalEvals = 0, localEvals = 0, refineLikTrees = 0, profileLikTrees = 0, ms_trees = 1,
		reportEveryEvals = 0, set_threads = 0, numGlobalSearches = 1;
size_t bestParsMapSize = 0;

double globalUpper = 5, globalLower = 1e-3, dataLnL = 0.0, bestGlobalSlLnL = -1000000.0, bestLocalSlLnL = -1000000.0, userLnL = 0.0, localSearchAbsTol = 1e-3;
bool skipGlobalSearch = false, bSFSmode = false, profileLikBool = false, abortNLopt = false, cmdLineInConfigFile = false, progressiveBounds = false,
		seedPRNGBool = false, nobSFSFile = false, printLikCorrFactor = false, startRandom = false, dataConvert = false, skipLocalSearch = false;
unsigned long int finalTableSize = 0, seedPRNG = 123456;

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

	if (poisTableSize) {
		free(poisTableScaleOrder);
		freed2int_matrix(blockLengthsMat, poisTableSize);
	}
}


//	Temporarily deactivated code for likelihood profiles (not to be confused with profile likelihoods!)
void profileLik(vector<double> MLEparVec) {

	double parLowerBound, parUpperBound, MLEparVal, MLEparLik;

	size_t parIdx = 0;
	//	Likelihood profiles for all parameters to be inferred (tbi)
	for (map<int, vector<int> >::iterator it = tbiMsCmdIdx.begin(); it != tbiMsCmdIdx.end(); it++) {
		stringstream tbiFileName;
		MLEparVal = MLEparVec[tbi2ParVec[it->first]];

/*
		double quarterRange = (upperBounds[parIdx] - lowerBounds[parIdx]) / 4;
		parUpperBound = MLEparVal + quarterRange;
		parLowerBound = MLEparVal - quarterRange;
		if (parLowerBound < lowerBounds[parIdx])
			parLowerBound = lowerBounds[parIdx];
		if (parUpperBound > upperBounds[parIdx])
			parUpperBound = upperBounds[parIdx];
*/
		parUpperBound = upperBounds[parIdx];
		parLowerBound = lowerBounds[parIdx];

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

		//	Likelihood profile intervals for this parameter
		vector<double> parRange;
		if (parLowerBound > 0)
			parRange = logspaced(parLowerBound, parUpperBound, 10);
		else if (parUpperBound < 0)
			parRange = negLogspaced(parLowerBound, parUpperBound, 10);
		else
			parRange = linspaced(parLowerBound, parUpperBound, 10);

		for (size_t j = 0; j < parRange.size(); j++) {
			if (parRange[j] == MLEparVal)
				outFile << scientific << MLEparVal << "\t" << MLEparLik << endl;
			else {
				for (size_t i = 0; i < it->second.size(); i++) {
					stringstream ststProfilePar;
					ststProfilePar << parRange[j];
					ststProfilePar >> ms_argv[it->second[i]];
				}
				double loglik = computeLik();
				if (ms_crash_flag)
					loglik = 2*dataLnL;
				outFile << scientific << parRange[j] << "\t" << loglik << endl;
			}
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


int getBrConfigNum(vector<int> brConfVec) {
	return intVec2BrConfig[brConfVec];
}


void process_tree_cond_bSFS (double ***onetreePoisTable) {
	//	TODO: need to look for alternatives to speed up mbSFS traversal
	for (int data = 0; data < mbSFSLen; data++) {
		for (size_t i = 0; i < dataConfigs[data].size(); i++) {
			for (int scale = 0; scale < dataScales[data]; scale++) {
				double jointPoisson = 1.0;

				for (int j = 0; j < brClass; j++)
					jointPoisson *= onetreePoisTable[data2PoisTableIdx[data][scale]][j][dataConfigs[data][i][j]];

				if (jointPoisson > (numeric_limits<double>::min()*ms_trees)) {
					selectConfigFreqs[data][i] += jointPoisson;
					trackSelectConfigsForInf[data][i] = 1;
				}
			}
		}
	}
}


void process_tree_exact_bSFS (double ***onetreePoisTable) {
	for (unsigned long int i = 0; i < finalTableSize; i++) {
	    double jointPoisson = 1.0;

	    vector<int> vec = getMutConfigVec(i);
		for (int j = 0; j < brClass; j++)
			jointPoisson *= onetreePoisTable[0][j][vec[j]];

		if (jointPoisson > (numeric_limits<double>::min()*ms_trees)) {
			allConfigs[i] = i;
			allConfigFreqs[i] += jointPoisson;
		}
	}
}


void calcBSFSTable() {

	ms_crash_flag = 0;
	{
		stringstream stst;
		stst << 1;
		stst >> ms_argv[2];
	}

	if (estimate == 1) {

/*
		int treesSampled = 0;
#pragma omp parallel for
		for (int trees = 0; trees < ms_trees; trees++) {
			if (!ms_crash_flag) {
				double **onetreePoisTable;
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
*/
	}
	else {

		if (estimate == 0) {
			allConfigs = vector<int>(finalTableSize,-1);
			allConfigFreqs = vector<double>(finalTableSize,0.0);
		}

		int sim_trees = ms_trees, crashLimit;
		crash_counter = sampledTrees = 0;
		if (ms_trees < 100)
			crashLimit = 1;
		else
			crashLimit = ms_trees/100;

		do {

#pragma omp parallel for shared(ms_crash_flag, crash_counter, sampledTrees)
			for (int trees = 0; trees < sim_trees; trees++) {
				if (!ms_crash_flag) {
					double ***onetreePoisTable;
					onetreePoisTable = d3matrix(poisTableSize,brClass, mutClass);
					// calling ms for sampling genealogies
					if (main_ms_ABLE(ms_argc, ms_argv, onetreePoisTable)) {
#pragma omp atomic
						++crash_counter;
					}
					else {
						if (bSFSmode || (estimate > 0))
							process_tree_cond_bSFS(onetreePoisTable);
						else if (estimate == 0)
							process_tree_exact_bSFS(onetreePoisTable);
#pragma omp atomic
						++sampledTrees;
					}

					if (crash_counter >= crashLimit)
						ms_crash_flag = 1;

					freed3matrix(onetreePoisTable, poisTableSize, brClass);
				}
			}

			if (sampledTrees == ms_trees)
				break;
			else if ((sampledTrees < ms_trees) && !ms_crash_flag)
				sim_trees = ms_trees - sampledTrees;

		} while(!ms_crash_flag);
	}
}


double computeLik() {

	double loglik = 0.0;

	if (bSFSmode || (estimate > 0)) {
		for (int data = 0; data < mbSFSLen; data++) {
			selectConfigFreqs.push_back(vector<double>(dataConfigFreqs[data].size(),0.0));
			trackSelectConfigsForInf.push_back(vector<int>(dataConfigFreqs[data].size(),0));
		}
	}

	//	calculating the bSFS config. probs.
	calcBSFSTable();

	if (!ms_crash_flag) {

		vector<double> loglikData(mbSFSLen, 0.0);
		for (int data = 0; data < mbSFSLen; data++) {

			vector<int> trackedConfigs (mbSFSLen, 0);
			if (bSFSmode || (estimate > 1))
				trackedConfigs[data] = accumulate(trackSelectConfigsForInf[data].begin(),trackSelectConfigsForInf[data].end(),0);
//			else if (estimate == 1)
//				trackedConfigs = trackSelectConfigs.size();

			if (bSFSmode) {
				if (nobSFSFile || profileLikBool) {
					for (size_t i = 0; i < dataConfigs[data].size(); i++) {
						if (selectConfigFreqs[data][i] != 0.0)
							loglikData[data] += log(selectConfigFreqs[data][i] / (ms_trees * dataScales[data])) * dataConfigFreqs[data][i];
					}
				}
				else {
					ofstream ofs(bSFSFile.c_str(),ios::out);

					for (size_t i = 0; i < dataConfigs[data].size(); i++) {
						if (selectConfigFreqs[data][i] != 0.0) {
							ofs << getMutConfigStr(dataConfigs[data][i]) << " : " << scientific << selectConfigFreqs[data][i] / (ms_trees * dataScales[data]) << endl;
							loglikData[data] += log(selectConfigFreqs[data][i] / (ms_trees * dataScales[data])) * dataConfigFreqs[data][i];
						}
					}

					ofs.close();
				}

				//	correct for the LnL if there are any data bSFS configs are unvisited
				if (trackedConfigs[data])
					loglikData[data] *= (double) dataConfigFreqs[data].size() / trackedConfigs[data];

				if (printLikCorrFactor)
					printf("Likelihood correction factor : %.6f\n", (double) dataConfigFreqs[data].size() / trackedConfigs[data]);
			}
			else if (estimate > 0) {
				for (size_t i = 0; i < dataConfigs[data].size(); i++) {
					if (selectConfigFreqs[data][i] != 0.0)
						loglikData[data] += log(selectConfigFreqs[data][i] / (ms_trees * dataScales[data])) * dataConfigFreqs[data][i];
				}

				//	correct for the LnL if there are any data bSFS configs are unvisited
				if (trackedConfigs[data])
					loglikData[data] *= (double) dataConfigFreqs[data].size() / trackedConfigs[data];
			}
		}

		if (!mbSFSLen && (estimate == 0)) {
			ofstream ofs(bSFSFile.c_str(),ios::out);
			for (size_t i = 0; i < allConfigs.size(); i++)
				if (allConfigs[i] >= 0)
					ofs << getMutConfigStr(allConfigs[i]) << " : " << scientific << allConfigFreqs[allConfigs[i]] / ms_trees << endl;
			ofs.close();

			allConfigs.clear();
			allConfigFreqs.clear();
		}

		loglik = accumulate(loglikData.begin(), loglikData.end(), 0.0);
	}

	trackSelectConfigsForInf.clear();
//	trackSelectConfigs.clear();
	selectConfigFreqs.clear();

	return loglik;
}


//	conversion from decimal to base-(maxPopSize+1)
void evalBranchConfigs() {

	int quo, rem, maxPopSize, totPopSum, count = 0, sumConfig;
	bool skipConfig;
	maxPopSize = totPopSum = sampledPops[0];

	if (kmax == 0)
		mutClass = 0;
	else
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


void readConfigFile() {

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
			if (tokens[0] == "ABLE") {
				cmdLineInConfigFile = true;
				cmdLine = tokens;
			}
			else if (tokens[0] == "pops") {
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
				else if (tokens[1] == "random") {
					startRandom = true;
				}
				else {
					stringstream stst(tokens[2]);
					stst >> val;
					tbiUserVal[atoi(tokens[1].substr(3).c_str())] = val;
				}
			}
			else if (tokens[0] == "datafile") {
				mbSFSLen = 0;
				for(unsigned int j = 1; j < tokens.size(); j++) {
					dataFile.push_back(tokens[j]);
					++mbSFSLen;
				}
			}
			else if (tokens[0] == "scale") {
				for(unsigned int j = 1; j < tokens.size(); j++) {
					dataScales.push_back(atoi(tokens[j].c_str()));
					uniqueScales.insert(atoi(tokens[j].c_str()));
				}
			}
			else if (tokens[0] == "datafile_format") {
				dataFileFormat = tokens[1];
			}
			else if (tokens[0] == "convert_data_to_bSFS") {
				dataConvert = true;
				data2bSFSFile = "data2bSFS.txt";
				if (tokens.size() > 1)
					data2bSFSFile = tokens[1];
			}
			else if (tokens[0] == "allele_type") {
				alleleType = tokens[1];
			}
			//	this option has been deprecated
			else if (tokens[0] == "estimate") {
				stringstream stst(tokens[1]);
				stst >> estimate;
			}
			else if (tokens[0] == "task") {
				if (tokens[1] == "exact_bSFS")
					estimate = 0;
				else if (tokens[1] == "conditional_bSFS")
					bSFSmode = true;
				else if (tokens[1] == "infer")
					estimate = 2;
			}
			else if (tokens[0] == "bSFS") {
				if (tokens.size() > 1)
					bSFSFile = tokens[1];
			}
			else if (tokens[0] == "kmax") {
				stringstream stst(tokens[1]);
				stst >> kmax;
			}
			else if (tokens[0] == "folded")
				foldBrClass = 1;
			else if (tokens[0] == "global_search")
				globalSearchAlg = tokens[1];
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
				skipGlobalSearch = true;
			}
			else if (tokens[0] == "skip_local_search") {
				skipLocalSearch = true;
			}
			else if (tokens[0] == "bounds") {
				int paramID = atoi(tokens[1].substr(3).c_str());
				for(unsigned int j = 2; j < tokens.size(); j++) {
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
			else if (tokens[0] == "retain_best_points") {
				if (tokens.size() > 1) {
					stringstream stst(tokens[1]);
					stst >> bestParsMapSize;
				}
				else
					bestParsMapSize = 10;
			}
			else if (tokens[0] == "global_search_iterations") {
				if (tokens.size() > 1) {
					stringstream stst(tokens[1]);
					stst >> numGlobalSearches;
				}
			}
			else if (tokens[0] == "global_bounds_refinement") {
				if (tokens.size() > 1) {
					if (tokens[1] == "progressive")
						progressiveBounds = true;
					else if (tokens[1] == "overall")
						progressiveBounds = false;
				}
			}
			else if (tokens[0] == "seed_PRNG") {
				stringstream stst(tokens[1]);
				stst >> seedPRNG;
				seedPRNGBool = true;
			}
			else if (tokens[0] == "refine_likelihoods") {
				if (tokens.size() > 1) {
					stringstream stst(tokens[1]);
					stst >> refineLikTrees;
				}
			}
			else if (tokens[0] == "profile_likelihoods") {
				bSFSmode = profileLikBool = true;
				if (tokens.size() > 1) {
					stringstream stst(tokens[1]);
					stst >> profileLikTrees;
				}
			}
			else if (tokens[0] == "report_likelihoods") {
				stringstream stst(tokens[1]);
				stst >> reportEveryEvals;
			}
			else if (tokens[0] == "start_likelihood") {
				stringstream stst(tokens[1]);
				stst >> userLnL;
			}
			else if (tokens[0] == "no_bSFS_file") {
				nobSFSFile = true;
			}
			else if (tokens[0] == "print_correction_factor") {
				printLikCorrFactor = true;
			}
			else if (tokens[0] == "set_ftol_abs") {
				stringstream stst(tokens[1]);
				stst >> localSearchAbsTol;
			}
			else if (tokens[0] == "set_threads") {
				stringstream stst(tokens[1]);
				stst >> set_threads;
			}
			else if (tokens[0] == "profile") {
				int paramID = atoi(tokens[1].substr(3).c_str());
				for(unsigned int j = 2; j < tokens.size(); j++) {
					double val;
					stringstream stst(tokens[j]);
					stst >> val;
					tbiProfilesGrid[paramID].push_back(val);
				}
			}
			else {
				cerr << "Unrecognised keyword \"" << tokens[0] << "\" found in the config file!" << endl;
				cerr << "Aborting ABLE..." << endl;
				exit(-1);
			}
		}
	}
	ifs.close();

}


void parseCmdLine(char* argv[]) {

	if (cmdLineInConfigFile)
		ms_argc = cmdLine.size();

	ms_argv = (char **)malloc( ms_argc*sizeof(char *) ) ;
	for(int i =0; i < ms_argc; i++)
		ms_argv[i] = (char *)malloc(30*sizeof(char) ) ;

	for (int i = 0; i < ms_argc; i++) {
		string param;
		if (cmdLineInConfigFile)
			param = cmdLine[i];
		else
			param = string(argv[i]);

		stringstream stst;
		if (param.substr(0,3) == "tbi") {
			if ((estimate == 2) || profileLikBool) {
				int paramID = atoi(param.substr(3).c_str());
				tbiMsCmdIdx[paramID].push_back(i);
				if (tbiUserVal.find(paramID) != tbiUserVal.end()) {
					stst << tbiUserVal[paramID];
					stst >> ms_argv[i];
				}
				else {
					//	check for tbi start values which are needed for a local search
					//	N.B. Likelihood profile code needs to re reviewed when activated in the future
					if (skipGlobalSearch && !startRandom) {
						cerr << "For this task you need to specify values for the \"tbi\" keywords using the \"start\" keyword" << endl;
						cerr << "Exiting ABLE...\n" << endl;
						free_objects();
						exit(-1);
					}

					double tmpPar = gsl_rng_uniform(PRNG);
					tbiUserVal[paramID] = tmpPar;
					stst << tmpPar;
					stst >> ms_argv[i];

					if (startRandom)
						setRandomPars[paramID] = true;
				}
			}
			else {
				cerr << "You cannot use \"tbi\" keywords for this task!" << endl;
				cerr << "Exiting ABLE...\n" << endl;
				free_objects();
				exit(-1);
			}
		}
		else {
			if (cmdLineInConfigFile)
				stst << cmdLine[i];
			else
				stst << argv[i];

			stst >> ms_argv[i];

			if (param == "-r") {
				stringstream length;
				if (cmdLineInConfigFile)
					length << cmdLine[i+2];
				else
					length << argv[i+2];

				length >> recLen;
			}
		}
	}

/*
	for(int i = 1; i < ms_argc; i++)
		printf("%s ",ms_argv[i]);
	printf("\n");
	free_objects();
	exit(-1);
*/
}


void checkConfigOptions() {
	if ((estimate == 2) || profileLikBool) {
		if (!tbiMsCmdIdx.size()) {
			cerr << "\nCannot proceed with inference" << endl;
			cerr << "\"tbi\" (To be Inferred) keywords are expected in the \"task infer\" mode" << endl;
			cerr << "Exiting ABLE...\n" << endl;
			free_objects();
			exit(-1);
		}
		else if (estimate == 2) {
			stringstream stst;
			if (globalTrees == 0)
				globalTrees = 10000;

			if (localTrees == 0)
				localTrees = 15000;
		}
	}

	if ((dataFile.size() == 1) && (dataScales.size() == 0)) {
		dataScales.push_back(1);
		uniqueScales.insert(1);
	}

	if (dataFile.size() != dataScales.size()) {
		cerr << "The respective \"scale\" of each specified dataset must be specified for using the mbSFS!" << endl;
		cerr << "Exiting ABLE...\n" << endl;
		free_objects();
		exit(-1);
	}
	else if ((estimate > 1) || bSFSmode) {
		//	enumerate and label the number of onetreePoisTables needed
		poisTableSize = accumulate(uniqueScales.begin(), uniqueScales.end(), 0);
		poisTableScaleOrder = (int *) malloc(poisTableSize * sizeof(int));
		//	calculate the block lengths for each "SUB"-onetreePoisTable of the mbSFS
		blockLengthsMat = d2int_matrix(poisTableSize, 3);
		int count = 0;
		for (set<int>::iterator it = uniqueScales.begin(); it != uniqueScales.end(); it++) {
			int prev = -1;
			for (int i = 0; i < *it; i++) {
				poisTableScaleOrder[count] = *it;

				if (recLen) {
					blockLengthsMat[count][0] = recLen / *it;
					blockLengthsMat[count][1] = prev + 1;
					blockLengthsMat[count][2] = prev + blockLengthsMat[count][0];
				}
				else {
					blockLengthsMat[count][0] = 0;
					blockLengthsMat[count][1] = 0;
					blockLengthsMat[count][2] = 0;
				}

				prev = blockLengthsMat[count][2];
				++count;
			}
		}

		//	data2PoisTableIdx : maps each dataset to its poisTable index/indices
		data2PoisTableIdx = vector<vector<int> > (dataScales.size(), vector <int>());
		for (size_t i = 0; i < dataScales.size(); i++)
			for (int j = 0; j < poisTableSize; j++)
				if (dataScales[i] == poisTableScaleOrder[j])
					data2PoisTableIdx[i].push_back(j);
	}
	else if (estimate == 0) {
		poisTableSize = 1;
		poisTableScaleOrder = (int *) malloc(poisTableSize * sizeof(int));
		poisTableScaleOrder[0] = 1;
		blockLengthsMat = d2int_matrix(poisTableSize, 3);
		if (recLen) {
			blockLengthsMat[0][0] = recLen;
			blockLengthsMat[0][1] = 0;
			blockLengthsMat[0][2] = recLen;
		}
		else {
			blockLengthsMat[0][0] = 0;
			blockLengthsMat[0][1] = 0;
			blockLengthsMat[0][2] = 0;
		}
	}

/*
	for (int block = 0; block < poisTableSize; block++)
		printf("Length: %d; start: %d; end: %d\n", blockLengthsMat[block][0], blockLengthsMat[block][1], blockLengthsMat[block][2]);
	free_objects();
	exit(-1);
*/
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

	if (((currState == GLOBAL) && (evalCount >= globalEvals))
			|| ((currState == LOCAL) && localEvals && (evalCount >= localEvals))) {
		AUGLAG.set_force_stop(2);
		return 0.0;
	}

	//	checking for simple non linear constraints between the free params
	for (map<int , int>::iterator it = parConstraints.begin(); it != parConstraints.end(); it++) {
		if (vars[tbi2ParVec[it->first]] >= vars[tbi2ParVec[it->second]]) {
			parConstraintPass = false;
			break;
		}
	}

	if (!parConstraintPass) {
/*
		printf("%5d ", evalCount);
		for (size_t i = 0; i < vars.size(); i++)
			printf("%.6f ", vars[i]);

		printf(" Trees: %d ", 0);
		printf(" Penalised LnL: %.6f\n", 10*dataLnL);
*/
		return 10*dataLnL;
	}

	//	constructing the ms command line
	for (map<int, vector<int> >::iterator it = tbiMsCmdIdx.begin(); it != tbiMsCmdIdx.end(); it++) {
		for (size_t i = 0; i < it->second.size(); i++) {
			stringstream stst;
			stst << vars[tbi2ParVec[it->first]];
			stst >> ms_argv[it->second[i]];
		}
	}
/*
	for(int i = 1; i < ms_argc; i++)
		printf("%s ",ms_argv[i]);
	printf("\n");
//		free_objects(); exit(-1);
*/

	//	MAIN COMPOSITE LIKELIHOOD COMPUTATION
	loglik = computeLik();

	if (ms_crash_flag) {
/*
		printf("%5d ", evalCount);
		for (size_t i = 0; i < vars.size(); i++)
			printf("%.6f ", vars[i]);

		printf(" Trees: %d ", 0);
		printf(" Penalised LnL: %.6f ", 10*dataLnL);
		printf(" ms CRASH!\n");
*/
		return 10*dataLnL;
	}

	//	pretty output
	++evalCount;
	printf("%6d ", evalCount);
	for (size_t i = 0; i < vars.size(); i++)
		printf("%.6f ", vars[i]);
//	printf(" Trees: %d ", ms_trees);
//	printf(" Sampled: %d ", sampledTrees);
	printf(" LnL: %.6f\n", loglik);

	if (loglik == 0.0) {
		abortNLopt = true;
		AUGLAG.set_force_stop(-5);
		return 0.0;
	}

	if ((currState == GLOBAL) && reportEveryEvals && !(evalCount % reportEveryEvals) && (evalCount < globalEvals)) {
		printf("\nReporting the best MLE after %d evaluations\n", evalCount);
		for (size_t i = 0; i < bestGlobalSPars.size(); i++)
			printf("%.6f ", bestGlobalSPars[i]);
		printf(" LnL: %.6f\n\n", bestGlobalSlLnL);
	}

	//	side stepping nlopt by storing the best LnL and parameters
	if (currState == GLOBAL) {
		if (loglik > bestGlobalSlLnL) {
			bestGlobalSlLnL = loglik;
			bestGlobalSPars = vars;
		}

		if (bestParsMapSize) {
			bestParsMap[loglik] = vars;
			if (bestParsMap.size() > bestParsMapSize)
				bestParsMap.erase(bestParsMap.begin()->first);
		}
	}
	else if ((currState == LOCAL) && (loglik > bestLocalSlLnL)) {
			bestLocalSlLnL = loglik;
			bestLocalSPars = vars;
	}

	return loglik;
}


double check_constraints(const vector<double> &vars, vector<double> &grad, void *data) {

	double consDiff = 1e-7;

	if (!grad.empty()) {
		cerr << "Cannot proceed with ABLE" << endl;
		cerr << "Gradient based optimization not yet implemented..." << endl;
		free_objects();
		exit(-1);
	}

	//	checking for simple non linear constraints between the free params
	for (map<int , int>::iterator it = parConstraints.begin(); it != parConstraints.end(); it++)
		consDiff += (vars[tbi2ParVec[it->first]] - vars[tbi2ParVec[it->second]]);

/*
	printf("%5d ", evalCount);
	for (size_t i = 0; i < vars.size(); i++)
		printf("%.6f ", vars[i]);

	printf(" diff: %f\n", consDiff);
*/

	return consDiff;
}


//	recursive exploration of user-specified parameter profiles
void exploreProfiles (size_t &varIdx) {
	for(unsigned int i = 0; i < tbiProfilesGrid[profileVarKey[varIdx]].size(); i++) {
		//	recursively constructing the parameter profile combination to be evaluated
		profileVars[varIdx] = tbiProfilesGrid[profileVarKey[varIdx]][i];
		if (varIdx+1 < profileVarKey.size())
			exploreProfiles(++varIdx);

		//	finished constructing the parameter profile combination to be evaluated
		else if (varIdx+1 == profileVarKey.size()) {
			//	constructing the ms command line
			for (map<int, vector<int> >::iterator it = tbiMsCmdIdx.begin(); it != tbiMsCmdIdx.end(); it++) {
				for (size_t i = 0; i < it->second.size(); i++) {
					stringstream stst;
					stst << profileVars[tbi2ParVec[it->first]];
					stst >> ms_argv[it->second[i]];
				}
			}

			bool parConstraintPass = true;
			//	checking for simple non linear constraints between the free params
			for (map<int , int>::iterator it = parConstraints.begin(); it != parConstraints.end(); it++) {
				if (profileVars[tbi2ParVec[it->first]] >= profileVars[tbi2ParVec[it->second]]) {
					parConstraintPass = false;
					break;
				}
			}

			//	pretty printing parameter profile combinations
			++evalCount;
			printf("%6d ", evalCount);
			for (size_t j = 0; j < profileVars.size(); j++)
				printf("%.6f ", profileVars[j]);

			//	IF parameter combination satisfies user-specified constraints
			if (parConstraintPass)
				printf(" LnL: %.6f\n", computeLik());
			else
				printf(" Point does not pass user-specified constraint...skipping!\n");
		}
	}
	--varIdx;
}


int main(int argc, char* argv[]) {

	string version = "0.1.x (Built on " + datestring + " at " + timestring + ")";

	cout << endl << endl;
	cout << "******************************************************************" << endl;
	cout << "*  This is ABLE version " << version << ".  *" << endl;
	cout << "*  ABLE is distributed under the CeCILL licence. See             *" << endl;
	cout << "*  http://www.cecill.info/index.en.html for more information.    *" << endl;
	cout << "*  © Champak Beeravolu Reddy 2015-now (champak.br@gmail.com)     *" << endl;
	cout << "******************************************************************" << endl;
	cout << endl << endl;

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

	readConfigFile();

	evalBranchConfigs();

	if ((estimate > 1) || bSFSmode || dataConvert) {
		if (dataFileFormat == "pseudo_MS") {
			readDataAsSeqBlocks("block_SNPs.txt", alleleType);

			//	convert data into bSFS and quit ABLE
			//	restricting this task to converting one dataset at a time
			if (dataConvert && (mbSFSLen == 1)) {
				ofstream ofs(data2bSFSFile.c_str(),ios::out);
				for (size_t i = 0; i < dataConfigs[0].size(); i++)
					ofs << getMutConfigStr(dataConfigs[0][i]) << " : " << setprecision(5) << scientific << dataConfigFreqs[0][i] << endl;
				ofs.close();

				cout << "Finished converting the data into a bSFS format...\n" << endl;
				return 0;
			}
		}
		else
			readDataAsbSFSConfigs();
	}

	parseCmdLine(argv);
	checkConfigOptions();

	int procs;
	if (set_threads > 0)
		procs = set_threads;
	else
		procs = omp_get_num_procs();

	omp_set_num_threads(procs);
	printf("Setting up %d threads...\n\n", procs);

	if (!seedPRNGBool)
		seedPRNG = hash(time(NULL), clock());

	for (int i = 0; i < procs; i++) {
		PRNGThreadVec.push_back(gsl_rng_alloc(gsl_rng_mt19937));
		gsl_rng_set(PRNGThreadVec[i], seedPRNG + i);
	}


//*********************************************************************************************************************************************
	//	i.e. task infer
	if (estimate == 2) {
		double maxLnL;
		bool endTaskInfer = false;

		vector<double> parVec;
		int count = 0;
		for (map<int, vector<int> >::iterator it = tbiMsCmdIdx.begin(); it != tbiMsCmdIdx.end(); it++) {
			//	Customising search bounds based on user specifications
			if ((tbiSearchBounds.find(it->first) != tbiSearchBounds.end()) && (tbiSearchBounds[it->first].size() >= 2)) {
				lowerBounds.push_back(tbiSearchBounds[it->first][0]);
				upperBounds.push_back(tbiSearchBounds[it->first][1]);
				//	if HARD upper/lower bounds have been specified
				if (tbiSearchBounds[it->first].size() == 4) {
					double hardLower = tbiSearchBounds[it->first][2],
						   hardUpper = tbiSearchBounds[it->first][3];
					if ((lowerBounds.back() < hardLower) || (hardUpper < upperBounds.back())) {
						cerr << "Please recheck your bounds for tbi" << it->first << ", the Upper/Lower hard bounds seem to have been misspecified!" << endl;
						cerr << "Aborting ABLE..." << endl;
						free_objects();
						exit(-1);
					}
					else {
						hardLowerBounds.push_back(hardLower);
						hardUpperBounds.push_back(hardUpper);
					}
				}
				else if (tbiSearchBounds[it->first].size() == 2) {
					hardLowerBounds.push_back(tbiSearchBounds[it->first][0]);
					hardUpperBounds.push_back(tbiSearchBounds[it->first][1]);
				}
			}
			else {
				lowerBounds.push_back(globalLower);
				upperBounds.push_back(globalUpper);
				hardLowerBounds.push_back(globalLower);
				hardUpperBounds.push_back(globalUpper);
			}
			tbi2ParVec[it->first] = count;
			ParVecToTbi[count] = it->first;

			//	if random user values have been specified
			double parVal;
			if (setRandomPars[it->first]) {
				//	breaks the PRNG dependence on the CPU clock in the case of a simultaneous start of jobs (e.g. on a cluster)
				//	and a random initialisation of the start values
				gsl_rng_set(PRNG, seedPRNG - 1);
				double tmpPar = gsl_rng_uniform(PRNG);
				parVal = tmpPar * (upperBounds[count] - lowerBounds[count]) + lowerBounds[count];
			}
			else
				parVal = tbiUserVal[it->first] * (upperBounds[count] - lowerBounds[count]) + lowerBounds[count];

			parVec.push_back(parVal);

			//	constructing the ms command line
			for (size_t i = 0; i < it->second.size(); i++) {
				stringstream stst;
				stst << parVec[tbi2ParVec[it->first]];
				stst >> ms_argv[it->second[i]];
			}
			++count;
		}

		//	Specifying the the Augmented Lagrangian algorithm (http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms#Augmented_Lagrangian_algorithm)
		void* data;
		AUGLAG = nlopt::opt(nlopt::AUGLAG, tbiMsCmdIdx.size());
		AUGLAG.set_lower_bounds(lowerBounds);
		AUGLAG.set_upper_bounds(upperBounds);
		AUGLAG.set_max_objective(optimize_wrapper_nlopt, NULL);
		AUGLAG.add_inequality_constraint(check_constraints, data, 1e-6);

		time(&likStartTime);

		//	If global search is not meant to be skipped... (more details the below here : http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms#Global_optimization)
		if (!skipGlobalSearch) {
			if (globalSearchAlg == "DIRECT") {
				opt = nlopt::opt(nlopt::GN_DIRECT, tbiMsCmdIdx.size());
				printf("Using the DIRECT algorithm for the global search...\n");
			}
			else if (globalSearchAlg == "CRS") {
				opt = nlopt::opt(nlopt::GN_CRS2_LM, tbiMsCmdIdx.size());
				opt.set_population(30*(tbiMsCmdIdx.size()+1));
				printf("Using the CONTROLLED RANDOM SEARCH WITH LOCAL MUTATION algorithm for the global search...\n");
			}
			else if (globalSearchAlg == "ISRES") {
				opt = nlopt::opt(nlopt::GN_ISRES, tbiMsCmdIdx.size());
				opt.set_population(30*(tbiMsCmdIdx.size()+1));
				printf("Using the IMPROVED STOCHASTIC RANKING EVOLUTION STRATEGY algorithm for the global search...\n");
			}
			else if (globalSearchAlg == "ESCH") {
				opt = nlopt::opt(nlopt::GN_ESCH, tbiMsCmdIdx.size());
				printf("Using Carlos Henrique da Silva Santos' EVOLUTIONARY algorithm for the global search...\n");
			}
			else {
				opt = nlopt::opt(nlopt::GN_DIRECT_NOSCAL, tbiMsCmdIdx.size());
				printf("Using the DIRECT_NOSCAL algorithm (i.e. without scaling) for the global search...\n");
			}

			int globalRecEvals = 2000 * tbiMsCmdIdx.size();
			if (globalEvals < globalRecEvals) {
				if (globalEvals)
					printf("\nThe specified number of GLOBAL search points with respect to the number of \"tbi\" parameters is below the recommended value.\nPlease consider increasing \"global_search_evals\"...\n");
				else
					globalEvals = globalRecEvals;
			}

			printf("\nStarting global search...\n");

			//	Global search with AUGLAG
			AUGLAG.set_local_optimizer(opt);

			//	Iterate over the specified number of global searches
			int globalIter = 1;
			while (globalIter <= numGlobalSearches) {
				ms_trees = globalTrees;
				evalCount = 0;
				abortNLopt = false;
				currState = GLOBAL;

				try {
					AUGLAG.optimize(parVec, maxLnL);
					throw abortNLopt;
				}
				catch (...) {
					if (abortNLopt) {
						cerr << "Something went wrong in the calculation of the LnL during the global search!" << endl;
						cerr << "Aborting ABLE..." << endl;

						for (size_t i = 0; i < parVec.size(); i++)
							printf("%.6f ", parVec[i]);
						printf(" LnL: %.6f\n\n", maxLnL);

						free_objects();
						exit(-1);
					}
				}

				if (bestParsMapSize && (globalIter < numGlobalSearches)) {
//					cout << "\nglobalIter: " << globalIter << endl;
//					cout << "bestLnL: " << bestParsMap.rbegin()->first << endl;
//					cout << "worstLnL: " << bestParsMap.begin()->first << endl;
//					for (map<double, vector<double> >::reverse_iterator it = bestParsMap.rbegin(); it != bestParsMap.rend(); it++) {
//						for (size_t i = 0; i < it->second.size(); i++)
//							printf("%.6f ", it->second[i]);
//						printf(" LnL: %.6f\n", it->first);
//					}
//					printf("\n");

					if (progressiveBounds) {
						printf("\nUsing the %d best global search results to automatically set bounds for the next GLOBAL search...\n", (int)bestParsMapSize);
						for (size_t i = 0; i < bestGlobalSPars.size(); i++) {
							set<double> sortParVals;
							for (map<double, vector<double> >::iterator it = bestParsMap.begin(); it != bestParsMap.end(); it++)
								sortParVals.insert(it->second[i]);

							double tmpLowerBound, tmpUpperBound;
							tmpLowerBound = *sortParVals.begin();
							tmpUpperBound = *sortParVals.rbegin();

//							printf("tbi%d: %.3f - %.3f; %.3f - %.3f; %.3f - %.3f \n", ParVecToTbi[i], tmpLowerBound, tmpUpperBound, lowerBounds[i], upperBounds[i], hardLowerBounds[i], hardUpperBounds[i]);

							//	IF tmpLowerBound is within 1% of the current lowerBounds[i]
							//	OR IF bestGlobalSPars[i] is within 5% of the current lowerBounds[i]
							//	then further DECREASE lowerBounds[i]
							if ((abs((tmpLowerBound - lowerBounds[i]) / lowerBounds[i]) < 0.01)
									|| (abs((bestGlobalSPars[i] - lowerBounds[i]) / lowerBounds[i]) < 0.05)) {
//								printf("Decreasing the lower bound for tbi%d...\n", ParVecToTbi[i]);
								if (lowerBounds[i] == hardLowerBounds[i])
									printf("WARNING. ABLE has attained the hard lower bound (%.3f) for tbi%d!\n", hardLowerBounds[i], ParVecToTbi[i]);
								lowerBounds[i] = (hardLowerBounds[i] + lowerBounds[i]) / 2;
							}
							//	else INCREASE the LOWER bound (using tmpLowerBound)
							else {
//								printf("Increasing the lower bound for tbi%d...\n", ParVecToTbi[i]);
								lowerBounds[i] = (tmpLowerBound + lowerBounds[i]) / 2;
							}


							//	IF tmpUpperBound is within 1% of the current upperBounds[i]
							//	OR IF bestGlobalSPars[i] is within 5% of the current upperBounds[i]
							//	then further INCREASE upperBounds[i]
							if ((abs((tmpUpperBound - upperBounds[i]) / upperBounds[i]) < 0.01)
									|| (abs((bestGlobalSPars[i] - upperBounds[i]) / upperBounds[i]) < 0.05)) {
//								printf("Increasing the upper bound for tbi%d...\n", ParVecToTbi[i]);
								if (upperBounds[i] == hardUpperBounds[i])
									printf("WARNING. ABLE has attained the hard upper bound (%.3f) for tbi%d!\n", hardUpperBounds[i], ParVecToTbi[i]);
								upperBounds[i] = (hardUpperBounds[i] + upperBounds[i]) / 2;
							}
							//	else DECREASE the UPPER bound (using tmpUpperBound)
							else {
//								printf("Decreasing the upper bound for tbi%d...\n", ParVecToTbi[i]);
								upperBounds[i] = (tmpUpperBound + upperBounds[i]) / 2;
							}
						}

						AUGLAG.set_lower_bounds(lowerBounds);
						AUGLAG.set_upper_bounds(upperBounds);

						int par = 0;
						for (map<int, vector<int> >::iterator it = tbiMsCmdIdx.begin(); it != tbiMsCmdIdx.end(); it++) {
							printf("tbi%d: %.3f - %.3f \n", it->first, lowerBounds[par], upperBounds[par]);
							++par;
						}
						printf("\n");
					}
					else
						printf("\nOnto the next GLOBAL search. The %d best global search results OVERALL will be retained.\n", (int)bestParsMapSize);

					parVec = bestGlobalSPars;
					maxLnL = bestGlobalSlLnL;
				}
				++globalIter;
			}

			if (skipLocalSearch) {
				parVec = bestGlobalSPars;
				maxLnL = bestGlobalSlLnL;
			}
		}

	    //	If local search is not meant to be skipped
		//	Specifying the the Subplex algorithm (http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms#Sbplx_.28based_on_Subplex.29)
		if (!skipLocalSearch) {
			if (!skipGlobalSearch) {
				if (numGlobalSearches > 1)
					printf("\nReporting the final global search MLE after %d iterations of %d evaluations each : \n", numGlobalSearches, evalCount);
				else
					printf("\nReporting the final global search MLE after %d evaluations : \n", evalCount);
				for (size_t i = 0; i < bestGlobalSPars.size(); i++)
					printf("%.6f ", bestGlobalSPars[i]);
				printf(" LnL: %.6f\n\n", bestGlobalSlLnL);

				if (localTrees < globalTrees) {
					cerr << "It is strongly advised to rerun the local search with a value of \"local trees\" >= \"global trees\"!" << endl;
					cerr << "Aborting ABLE..." << endl;
					free_objects();
					exit(-1);
				}

				parVec = bestGlobalSPars;
				bestLocalSlLnL = bestGlobalSlLnL;
				printf("\nUsing the global search MLE as the starting point for a refined local search...\n\n");

				if (bestParsMapSize) {
//					cout << "bestLnL: " << bestParsMap.rbegin()->first << endl;
//					cout << "worstLnL: " << bestParsMap.begin()->first << endl;
//					for (map<double, vector<double> >::reverse_iterator it = bestParsMap.rbegin(); it != bestParsMap.rend(); it++) {
//						for (size_t i = 0; i < it->second.size(); i++)
//							printf("%.6f ", it->second[i]);
//						printf(" LnL: %.6f\n", it->first);
//					}
//					printf("\n");

					printf("Using the %d best global search results to automatically set bounds for the local search...\n", (int)bestParsMapSize);
					for (size_t i = 0; i < bestGlobalSPars.size(); i++) {
						set<double> sortParVals;
						for (map<double, vector<double> >::iterator it = bestParsMap.begin(); it != bestParsMap.end(); it++)
							sortParVals.insert(it->second[i]);

						double tmpLowerBound, tmpUpperBound;
						tmpLowerBound = *sortParVals.begin();
						tmpUpperBound = *sortParVals.rbegin();

//						printf("tbi%d: %.3f - %.3f; %.3f - %.3f; %.3f - %.3f \n", ParVecToTbi[i], tmpLowerBound, tmpUpperBound, lowerBounds[i], upperBounds[i], hardLowerBounds[i], hardUpperBounds[i]);

						//	IF tmpLowerBound is within 1% of the current lowerBounds[i]
						//	OR IF bestGlobalSPars[i] is within 5% of the current lowerBounds[i]
						//	then further DECREASE lowerBounds[i]
						if ((abs((tmpLowerBound - lowerBounds[i]) / lowerBounds[i]) < 0.01)
								|| (abs((bestGlobalSPars[i] - lowerBounds[i]) / lowerBounds[i]) < 0.05)) {
//							printf("Decreasing the lower bound for tbi%d...\n", ParVecToTbi[i]);
							if (lowerBounds[i] == hardLowerBounds[i])
								printf("WARNING. ABLE has attained the hard lower bound (%.3f) for tbi%d!\n", hardLowerBounds[i], ParVecToTbi[i]);
							lowerBounds[i] = (hardLowerBounds[i] + lowerBounds[i]) / 2;
						}
						//	else INCREASE the LOWER bound (using tmpLowerBound)
						else {
//							printf("Increasing the lower bound for tbi%d...\n", ParVecToTbi[i]);
							lowerBounds[i] = (tmpLowerBound + lowerBounds[i]) / 2;
						}


						//	IF tmpUpperBound is within 1% of the current upperBounds[i]
						//	OR IF bestGlobalSPars[i] is within 5% of the current upperBounds[i]
						//	then further INCREASE upperBounds[i]
						if ((abs((tmpUpperBound - upperBounds[i]) / upperBounds[i]) < 0.01)
								|| (abs((bestGlobalSPars[i] - upperBounds[i]) / upperBounds[i]) < 0.05)) {
//							printf("Increasing the upper bound for tbi%d...\n", ParVecToTbi[i]);
							if (upperBounds[i] == hardUpperBounds[i])
								printf("WARNING. ABLE has attained the hard upper bound (%.3f) for tbi%d!\n", hardUpperBounds[i], ParVecToTbi[i]);
							upperBounds[i] = (hardUpperBounds[i] + upperBounds[i]) / 2;
						}
						//	else DECREASE the UPPER bound (using tmpUpperBound)
						else {
//							printf("Decreasing the upper bound for tbi%d...\n", ParVecToTbi[i]);
							upperBounds[i] = (tmpUpperBound + upperBounds[i]) / 2;
						}
					}

					AUGLAG.set_lower_bounds(lowerBounds);
					AUGLAG.set_upper_bounds(upperBounds);

					int par = 0;
					for (map<int, vector<int> >::iterator it = tbiMsCmdIdx.begin(); it != tbiMsCmdIdx.end(); it++) {
						printf("tbi%d: %.3f - %.3f \n", it->first, lowerBounds[par], upperBounds[par]);
						++par;
					}
					printf("\n");
				}
			}
			//	Local search without a prior global search
			else {
				printf("Skipping global search!\n");
				printf("\nUsing the user-specified/default values as a starting point for a local search...\n\n");
			}

			int localRecEvals = 1000 * tbiMsCmdIdx.size();
			if (localEvals < localRecEvals) {
				if (localEvals)
					printf("\nThe specified number of LOCAL search points with respect to the number of \"tbi\" parameters is below the recommended value.\nPlease consider increasing \"local_search_evals\"...\n");
				else
					localEvals = localRecEvals;
			}

			printf("\nStarting local search...\n");

			//	setting up local search params
			ms_trees = localTrees;
			evalCount = 0;
			abortNLopt = false;
			currState = LOCAL;

			vector<double> localSearchPerturb;
			for (size_t param = 0; param < lowerBounds.size(); param++)
				localSearchPerturb.push_back((upperBounds[param]-lowerBounds[param])/4);

			local_opt = nlopt::opt(nlopt::LN_SBPLX, tbiMsCmdIdx.size());
//			local_opt.set_xtol_rel(1e-2);
			local_opt.set_ftol_abs(localSearchAbsTol);
			local_opt.set_initial_step(localSearchPerturb);

			//	Local search with AUGLAG
			AUGLAG.set_local_optimizer(local_opt);

			try {
				AUGLAG.optimize(parVec, maxLnL);
				throw abortNLopt;
			}
			catch (...) {
				if (abortNLopt) {
					cerr << "Something went wrong in the calculation of the LnL during the local search!" << endl;
					cerr << "Aborting ABLE..." << endl;

					for (size_t i = 0; i < parVec.size(); i++)
						printf("%.6f ", parVec[i]);
					printf(" LnL: %.6f\n\n", maxLnL);

					free_objects();
					exit(-1);
				}
			}

			//	IF the local search started with global search results
			if (!skipGlobalSearch) {
				if (bestLocalSlLnL <= bestGlobalSlLnL) {
					printf("\nIgnoring local search results as they did not improve over the global search MLE...\n");
					parVec = bestGlobalSPars;
					maxLnL = bestGlobalSlLnL;
				}
				else {
					printf("\nFound a better local search MLE after %d evaluations\n", evalCount);
					parVec = bestLocalSPars;
					maxLnL = bestLocalSlLnL;
				}
			}
			//	IF user specified a cutoff LnL and/or no global search
			else {
				if ((userLnL != 0.0) && (maxLnL <= userLnL)) {
					printf("\nIgnoring local search results as they did not improve on the user-specified LnL after %d evaluations\n", evalCount);
					endTaskInfer = true;
				}
				else {
					printf("\nFound the local search MLE after %d evaluations\n", evalCount);
					parVec = bestLocalSPars;
					maxLnL = bestLocalSlLnL;
				}
			}
		}

		//	If ABLE didn't terminate prematurely end (i.e. if user specified LnL cutoff was not met during local search)
		if (!endTaskInfer) {
			//	Refining MLE after global/local search
			if ((!skipGlobalSearch || !skipLocalSearch) && refineLikTrees) {
				printf("\nRefining the likelihood at the MLE using %d genealogies...\n", refineLikTrees);
				for (map<int, vector<int> >::iterator it = tbiMsCmdIdx.begin(); it != tbiMsCmdIdx.end(); it++) {
					for (size_t i = 0; i < it->second.size(); i++) {
						stringstream stst;
						stst << parVec[tbi2ParVec[it->first]];
						stst >> ms_argv[it->second[i]];
					}
				}

				currState = OTHER;
				ms_trees = refineLikTrees;
				maxLnL = computeLik();
			}

			//	Final MLE output after global/local search and followed by LnL refinement (if specified)
			printf("\n\nFound a maximum at ");
			for (size_t i = 0; i < parVec.size(); i++)
				printf("%.6f ", parVec[i]);
			printf(" LnL: %.6f\n\n", maxLnL);
		}

		time(&likEndTime);
		printf("\nTime taken for computation : %.0f seconds\n\n", float(likEndTime - likStartTime));
	}

//*********************************************************************************************************************************************
	//	i.e. task conditional_bSFS or for profile_likelihoods
	else if ((estimate == 1) || bSFSmode) {
		currState = OTHER;
		double loglik;

		if (bSFSmode) {
			time(&likStartTime);
			bool zeroTrees = false;

			if (profileLikBool) {
				if (!profileLikTrees) {
					cerr << "\nPlease specify a non-zero number of trees/genealogies (to be used for computing the profile likelihoods) on the command line!" << endl;
					zeroTrees = true;
				}
				else if (tbiMsCmdIdx.size() != tbiProfilesGrid.size())
					zeroTrees = true; // using it as a flag for a premature exit
				else {
					ms_trees = profileLikTrees;

					int count = 0;
					profileVars = vector<double>(tbiMsCmdIdx.size(),0.0);
					for (map<int, vector<int> >::iterator it = tbiMsCmdIdx.begin(); it != tbiMsCmdIdx.end(); it++) {
						tbi2ParVec[it->first] = count;
						profileVarKey.push_back(it->first);
						++count;
					}

					evalCount = 0;
					size_t varIdx = 0;
					printf("\nEvaluating profile likelihood combinations...\n");
					exploreProfiles(varIdx);
				}
			}
			else {
				if (!ms_trees) {
					cerr << "\nPlease specify non-zero number of trees/genealogies (to be used for computing the conditional_bSFS) on the command line!" << endl;
					zeroTrees = true;
				}
				else {
					printf("Evaluating point likelihood at : \n");
					for (int i = 0; i < ms_argc; i++)
						cout << ms_argv[i] << " ";
					cout << endl;

					{
						stringstream stst;
						stst << ms_argv[2];
						stst >> ms_trees;
					}

					loglik = computeLik();
				}
			}

			time(&likEndTime);

			if (zeroTrees) {
				cerr << "Exiting ABLE...\n" << endl;
				free_objects();
				exit(-1);
			}
			else {
				if (ms_crash_flag) {
					cerr << "\nABLE failed to simulate genealogies with the demographic parameters that have been specified!" << endl;
					cerr << "Please consider changing them or contact the author Champak B. Reddy (champak.br@gmail.com)" << endl;
					cerr << "Exiting ABLE...\n" << endl;

					free_objects();
					exit(-1);
				}
				else if (!profileLikBool)
					printf(" LnL: %.6f (Trees sampled : %d)\n", loglik, sampledTrees);

				printf("\nTime taken for computation : %.0f seconds\n\n", float(likEndTime - likStartTime));
			}
		}
//		temporary deactivation : need to review code in calcBSFSTable for this option
//		else {
//			testLik.open("logliks.txt",ios::out);
//			testConfig.open("propConfigs.txt",ios::out);
//			time(&likStartTime);
//			loglik = computeLik();
//			time(&likEndTime);
//			testLik.close();
//			testConfig.close();
//		}
	}

//*********************************************************************************************************************************************
	//	i.e. task exact_bSFS
	else if (estimate == 0) {

		if (!mutClass) {
			cerr << "\nThe \"kmax\" option needs to be specified in the config file for calculating the expected bSFS!" << endl;
			cerr << "Exiting ABLE...\n" << endl;
			free_objects();
			exit(-1);
		}

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
			printf("\nTime taken for computation : %.0f seconds\n", float(likEndTime - likStartTime));
		}
	}

	free_objects();

	return 0;
}




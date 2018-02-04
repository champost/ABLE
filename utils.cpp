/***************************************************************************
© Champak Beeravolu Reddy 2016-now

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
 * utils.cpp
 *
 *  Created on: 4 Jan 2016
 *      Author: champost
 */

#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <map>
#include <fstream>
#include <iostream>
#include <omp.h>

#include "utils.h"

//	Linearly spaced points between [a,b]
vector<double> linspaced(double a, double b, int n) {
    vector<double> array;
    if ((n == 0) || (n == 1) || (a == b))
    	array.push_back(b);
    else if (n > 1) {
		double step = (b - a) / (n - 1);
		int count = 0;
		while(count < n) {
			array.push_back(a + count*step);
			++count;
		}
    }
    return array;
}


//	Log-spaced points between [a,b] ONLY for (0 < a,b)
vector<double> logspaced(double a, double b, int n) {
    vector<double> array;
    if ((a > 0) && (b > 0)) {
        if ((n == 0) || (n == 1) || (a == b))
        	array.push_back(b);
        else if (n > 1) {
            double step = pow(b/a, 1.0/(n-1));
            int count = 0;
        	while(count < n) {
        		array.push_back(a * pow(step, count));
        		++count;
        	}
        }
    }
    return array;
}


//	Log-spaced points between [a,b] ONLY for (a,b < 0)
vector<double> negLogspaced(double a, double b, int n) {
    vector<double> array;
    double posA = -a, posB = -b;
    if ((posB > 0) && (posA > 0)) {
        if ((n == 0) || (n == 1) || (posB == posA))
        	array.push_back(b);
        else if (n > 1) {
            double step = pow(posA/posB, 1.0/(n-1));
            int count = 0;
        	while(count < n) {
        		array.push_back(-posB * pow(step, count));
        		++count;
        	}
        }
    }
    reverse(array.begin(),array.end());
    return array;
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


void readDataAsSeqBlocks(string alleleType, bool outSNPs) {

	if (procs > 1) {
		outSNPs = false;
		printf("WARNING. No \"block_SNPs.txt\" will be created when using multiple CPU cores for converting data into the bSFS!\n");
	}

	int dataKmax = 0;
	dataConfigs = vector<vector<vector<int> > > (mbSFSLen, vector<vector<int> >());
	dataConfigFreqs = vector<vector<double> > (mbSFSLen, vector<double>());
	string outSNPsFile = "block_SNPs.txt";

	ofstream ofs;
	if (mbSFSLen == 1 && outSNPs)
		ofs.open(outSNPsFile.c_str(),ios::out);

	for (int data = 0; data < mbSFSLen; data++) {
		string line;
		ifstream ifs(dataFile[data].c_str(),ios::in);
		int nblocks = 0, configKmax;

		vector< vector<string> > blockMat;
		while (getline(ifs,line) && !ifs.eof()) {
			if (line[0] == '/') {
				blockMat.push_back(vector<string> ());

				//	skip header until the start of the block
				bool startBlock = false;
				while (getline(ifs,line) &&  !ifs.eof()) {
					if ((alleleType == "genotype") && ((line[0] == 'A') or (line[0] == 'T') or (line[0] == 'G') or (line[0] == 'C') or (line[0] == 'N'))) {
						startBlock = true;
						blockMat[nblocks].push_back(line);
					}
					else if ((alleleType == "binary") && ((line[0] == '0') or (line[0] == '1') or (line[0] == 'N'))) {
						startBlock = true;
						blockMat[nblocks].push_back(line);
					}
					else if (startBlock or (line.length() == 0))
						break;
				}
				++nblocks;
			}
		}
		ifs.close();

		vector < map<vector<int>, int> > finalTableMap = vector < map<vector<int>, int> > (procs, map<vector<int>, int> ());
#pragma omp parallel for
		for (size_t block = 0; block < blockMat.size(); block++) {
			int blockSize = 0, segSites = 0;
			bool monoMorphicBlock = true;
			vector<int> mutConfigVec(allBrClasses,0), foldedMutConfigVec(brClass,0);

			if (blockMat[block].size())
				blockSize = blockMat[block][0].size();

			for (int nuc = 0; nuc < blockSize; nuc++) {

				bool maskChar = false;
				vector<int> segCountVec;
				map<int, map<char, int> > popwiseAlleleCounts;
				map<char, bool> isAllele;
				int sample = 0;
				for (size_t pop = 0; pop < sampledPops.size(); pop++) {
					for (int popSam = 0; popSam < sampledPops[pop]; popSam++) {
						if (blockMat[block][sample][nuc] == 'N') {
							maskChar = true;
							pop = sampledPops.size() + 1;
							break;
						}
						++popwiseAlleleCounts[pop][blockMat[block][sample][nuc]];
						isAllele[blockMat[block][sample][nuc]] = true;
						++sample;
					}
				}
				if (!maskChar) {
					//	testing for bi-allelic SNP's
					int alleleCount = 0;
					char chooseAnAllele;
					for (map<char, bool>::iterator it = isAllele.begin(); it != isAllele.end(); it++) {
						if (it->second) {
							++alleleCount;
							if (alleleType == "genotype")
								chooseAnAllele = it->first;
							else if ((alleleType == "binary") && it->first == '1')
								chooseAnAllele = it->first;
						}
					}

					if (alleleCount == 2) {
						for (size_t pop = 0; pop < sampledPops.size(); pop++)
							segCountVec.push_back(popwiseAlleleCounts[pop][chooseAnAllele]);
						++mutConfigVec[getBrConfigNum(segCountVec)];

						//	NB. Tri/Quadri-allelic nucleotide positions are treated as monomorphic
						monoMorphicBlock = false;
					}

					if (alleleCount > 1)
						++segSites;
				}
			}
			if (mbSFSLen == 1 && outSNPs)
				ofs << segSites << endl;

			if (!monoMorphicBlock) {
				//	folding the multi-dimensional SFS
				for(int thisClass = 1; thisClass <= brClass; thisClass++) {
					foldedMutConfigVec[thisClass-1] = mutConfigVec[thisClass-1] + (foldBrClass * mutConfigVec[allBrClasses-thisClass]);
					if (foldBrClass && ((thisClass-1) == (allBrClasses-thisClass)))
						foldedMutConfigVec[thisClass-1] /= 2;
				}

				//	taking care of the Kmax
				for(int branch = 0; branch < brClass; branch++) {
					if ((mutClass > 0) && (foldedMutConfigVec[branch] > (mutClass - 2)))
						foldedMutConfigVec[branch] = mutClass - 1;
				}
			}
			++finalTableMap[omp_get_thread_num()][foldedMutConfigVec];
		}

		if (mbSFSLen == 1 && outSNPs)
			ofs.close();


		for (int thread = 1; thread < procs; thread++)
			for (map<vector<int>, int>::iterator it = finalTableMap[thread].begin(); it != finalTableMap[thread].end(); it++)
				finalTableMap[0][it->first] += it->second;

		for (map<vector<int>, int>::iterator it = finalTableMap[0].begin(); it != finalTableMap[0].end(); it++) {

	//		printf("%s : %.5e\n", getMutConfigStr(it->first).c_str(), (double) it->second/nblocks);

			dataConfigs[data].push_back(it->first);
			dataConfigFreqs[data].push_back((double) it->second/nblocks);

			configKmax = *max_element(it->first.begin(),it->first.end());
			if (configKmax > dataKmax)
				dataKmax = configKmax;
		}
	}

	dataLnL = 0.0;
	for (int data = 0; data < mbSFSLen; data++) {
		for (size_t i = 0; i < dataConfigs[data].size(); i++)
			dataLnL += log(dataConfigFreqs[data][i]) * dataConfigFreqs[data][i];
	}
	bestGlobalSlLnL = bestLocalSlLnL = 100000*dataLnL;

	if (!mutClass)
		mutClass = dataKmax;
}


void readDataAsbSFSConfigs() {

	int dataKmax = 0;
	dataConfigs = vector<vector<vector<int> > > (mbSFSLen, vector<vector<int> >());
	dataConfigFreqs = vector<vector<double> > (mbSFSLen, vector<double>());

	for (int data = 0; data < mbSFSLen; data++) {
		string line, del, keyVec;
		vector<string> tokens;
		vector<int> config;
		double val;
		int configKmax;

		ifstream ifs(dataFile[data].c_str(),ios::in);
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
			dataConfigs[data].push_back(config);
			dataConfigFreqs[data].push_back(val);

			configKmax = *max_element(config.begin(),config.end());
			if (configKmax > dataKmax)
				dataKmax = configKmax;

			config.clear();
		}
		ifs.close();
	}

	dataLnL = 0.0;
	for (int data = 0; data < mbSFSLen; data++) {
		for (size_t i = 0; i < dataConfigs[data].size(); i++)
			dataLnL += log(dataConfigFreqs[data][i]) * dataConfigFreqs[data][i];
	}
	bestGlobalSlLnL = bestLocalSlLnL = 100000*dataLnL;

	if (!mutClass)
		mutClass = dataKmax;
}




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
 * utils.h
 *
 *  Created on: 4 Jan 2016
 *      Author: champost
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <climits> // NOT the C++ version of limits... but defines UCHAR_MAX...

using namespace std;

typedef unsigned long uint32;

vector<double> linspaced(double a, double b, int n);
vector<double> logspaced(double a, double b, int n);
vector<double> negLogspaced(double a, double b, int n);
void Tokenize(const string& str, vector<string>& tokens, const string& delimiters);
void TrimSpaces(string& str);
void readDataAsSeqBlocks(string outSNPsFile, string alleleType);
void readDataAsbSFSConfigs();

extern vector<string> dataFile;
extern int brClass, mutClass, foldBrClass, allBrClasses, mbSFSLen;
extern double dataLnL, bestGlobalSlLnL, bestLocalSlLnL;
extern vector<int> sampledPops;
extern vector<vector<double> > dataConfigFreqs;
extern vector<vector<vector<int> > > dataConfigs;

extern int getBrConfigNum(vector<int> brConfVec);
extern string getMutConfigStr(vector<int> configVec);


inline uint32 hash( time_t t, clock_t c )
{
	// Get a uint32 from t and c
	// Better than uint32(x) in case x is floating point in [0,1]
	// Based on code by Lawrence Kirby (fred@genesis.demon.co.uk)

	static uint32 differ = 0;  // guarantee time-based seeds will change

	uint32 h1 = 0;
	unsigned char *p = (unsigned char *) &t;
	for( size_t i = 0; i < sizeof(t); ++i )
	{
		h1 *= UCHAR_MAX + 2U;
		h1 += p[i];
	}
	uint32 h2 = 0;
	p = (unsigned char *) &c;
	for( size_t j = 0; j < sizeof(c); ++j )
	{
		h2 *= UCHAR_MAX + 2U;
		h2 += p[j];
	}
	return ( h1 + differ++ ) ^ h2;
}

#endif /* UTILS_H_ */

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



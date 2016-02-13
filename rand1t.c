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


/***************************************************************************

This file was downloaded from https://webshare.uchicago.edu/users/rhudson1/Public/ms.folder/ms.tar.gz on Jan 14 2015
and thereafter modified by Champak Beeravolu Reddy (champak.br@gmail.com)

***************************************************************************/


/*  Link in this file for random number generation using drand48() */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_rng.h>

extern gsl_rng * prng;

         double
ran1()
{
/*
        double drand48();
        return( drand48() );

        int i;
    	for (i = 0; i < 10; i++) {
    		printf("%.5e\t%.5e\n", drand48(), ran2());
    	}
    	exit(-1);
*/
        	 return gsl_rng_uniform(prng);
/*
             int rand();
             return rand()/(RAND_MAX+1.0);
*/
}


/*
	void seedit( char *flag )
{
	FILE *fopen();
	unsigned short seedv[3], *seed48();

	if( flag[0] == 's' ){
	  seedv[0] =  (unsigned short) time(NULL) ;
            seedv[1] = 27011; seedv[2] = 59243;
          seed48( seedv );

       printf("\n%d %d %d\n", seedv[0], seedv[1], seedv[2] );
	}
}


	int
commandlineseed( char **seeds)
{
	unsigned short seedv[3], *seed48();

	seedv[0] = atoi( seeds[0] );
	seedv[1] = atoi( seeds[1] );
	seedv[2] = atoi( seeds[2] );
//	printf("\n%d %d %d\n", seedv[0], seedv[1], seedv[2] );

	seed48(seedv);
	return(3);
}
*/


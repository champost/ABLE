/*  Link in this file for random number generation using drand48() */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

         double
ran1()
{
//        double drand48();
//        return( drand48() );

//        int i;
//    	for (i = 0; i < 10; i++) {
//    		printf("%.5e\t%.5e\n", drand48(), ran2());
//    	}
//    	exit(-1);
        	 double ranMT();
        	 return ranMT();
}


//	void seedit( char *flag )
//{
//	FILE *fopen();
//	unsigned short seedv[3], *seed48();
//
//	if( flag[0] == 's' ){
//	  seedv[0] =  (unsigned short) time(NULL) ;
//            seedv[1] = 27011; seedv[2] = 59243;
//          seed48( seedv );
//
//       printf("\n%d %d %d\n", seedv[0], seedv[1], seedv[2] );
//	}
//}


//	int
//commandlineseed( char **seeds)
//{
//	unsigned short seedv[3], *seed48();
//
//	seedv[0] = atoi( seeds[0] );
//	seedv[1] = atoi( seeds[1] );
//	seedv[2] = atoi( seeds[2] );
////	printf("\n%d %d %d\n", seedv[0], seedv[1], seedv[2] );
//
//	seed48(seedv);
//	return(3);
//}


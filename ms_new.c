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


/***** ms.c     ************************************************
*
*       Generates samples of gametes ( theta given or fixed number 
*						of segregating sites.)
*	Usage is shown by typing ms without arguments.   
        usage: ms nsam howmany  -t  theta  [options]
		or
	       ms nsam howmany -s segsites  [options] 

	   nsam is the number of gametes per sample.
	   howmany is the number of samples to produce.
	   With -t the numbers of segregating sites will randomly vary 
		from one sample to the next.
	   with -s segsites,  the number of segregating sites will be
		segsites in each sample.

           Other options: See msdoc.pdf or after downloading and compiling, type ms<CR>.


*	  Arguments of the options are explained here:

           npop:  Number of subpopulations which make up the total population
           ni:  the sample size from the i th subpopulation (all must be 
		specified.) The output will have the gametes in order such that
		the first n1 gametes are from the first island, the next n2 are
		from the second island, etc.
           nsites: number of sites between which recombination can occur.
           theta: 4No times the neutral mutation rate 
           rho: recombination rate between ends of segment times 4No
	   f: ratio of conversion rate to recombination rate. (Wiuf and Hein model.)
	   track_len:  mean length of conversion track in units of sites.  The 
		       total number of sites is nsites, specified with the -r option.
           mig_rate: migration rate: the fraction of each subpop made up of
                 migrants times 4No. 
           howmany: howmany samples to generate.

	Note:  In the above definition, No is the total diploid population if
		npop is one, otherwise, No is the diploid population size of each
		subpopulation. 
	A seed file called "seedms" will be created  if it doesn't exist. The
		seed(s) in this file will be modified by the program. 
		So subsequent runs
		will produce new output.  The initial contents of seedms will be
		printed on the second line of the output.
        Output consists of one line with the command line arguments and one
	 	line with the seed(s).
		The samples appear sequentially following that line.
		Each sample begins with "//", then the number of segregating sites, the positions
		of the segregating sites (on a scale of 0.0 - 1.0). On the following
		lines are the sampled gametes, with mutants alleles represented as
		ones and ancestral alleles as zeros.
	To compile:  cc -o ms  ms.c  streec.c  rand1.c -lm
		or:  cc -o ms ms.c streec.c rand2.c -lm
	 (Of course, gcc would be used instead of cc on some machines.  And -O3 or 
		some other optimization switches might be usefully employed with some 
		compilers.) ( rand1.c uses drand48(), whereas rand2.c uses rand() ).

*
*   Modifications made to combine ms and mss on 25 Feb 2001
*	Modifications to command line options to use switches  25 Feb 2001
*	Modifications to add // before each sample  25 Feb 2001
	Modifications to add gene conversion 5 Mar 2001
	Added demographic options -d  13 Mar 2001
	Changed ran1() to use rand(). Changed seed i/o to accomodate this change. 20 April.
	Changed cleftr() to check for zero rand() .13 June 2001
	Move seed stuff to subroutine seedit()  11 July 2001
	Modified streec.c to handle zero length demographic intervals 9 Aug 2001
	Corrected problem with negative growth rates (Thanks to D. Posada and C. Wiuf) 13 May 2002
	Changed sample_stats.c to output thetah - pi rather than pi - thetah.  March 8 2003.
	Changed many command line options, allowing arbitrary migration matrix, and subpopulation
	   sizes.  Also allows parameters to come from a file. Option to output trees.  Option to
	   split and join subpopulations.   March 8, 2003. (Old versions saved in msold.tar ).
	!!! Fixed bug in -en option.  Earlier versions may have produced garbage when -en ... used. 9 Dec 2003
	Fixed bug which resulted in incorrect results for the case where
             rho = 0.0 and gene conversion rate > 0.0. This case was not handled
	    correctly in early versions of the program. 5 Apr 2004.  (Thanks to
	    Vincent Plagnol for pointing out this problem.) 
	Fixed bug in prtree().  Earlier versions may have produced garbage when the -T option was used.
		 1 Jul 2004.
	Fixed bug in -e. options that caused problems with -f option  13 Aug 2004.
	Fixed bug in -es option, which was a problem when used with -eG. (Thanks again to V. Plagnol.) 6 Nov. 2004
	Added -F option:  -F minfreq  produces output with sites with minor allele freq < minfreq filtered out.  11 Nov. 2004.
	Fixed bug in streec.c (isseg() ).  Bug caused segmentation fault, crash on some machines. (Thanks
	    to Melissa Jane Todd-Hubisz for finding and reporting this bug.)
	Added -seeds option 4 Nov 2006
	Added "tbs" arguments feature 4 Nov 2006
	Added -L option.  10 May 2007
	Changed -ej option to set Mki = 0 pastward of the ej event.  See msdoc.pdf.  May 19 2007.
	fixed bug with memory allocation when using -I option. This caused problems expecially on Windows
          machines.  Thanks to several people, including Vitor Sousa and Stephane De Mita for help on this one.
          Oct. 17, 2007.
     Modified pickb() and pickbmf() to eliminate rare occurrence of fixed mutations Thanks to J. Liechty and K. Thornton. 4 Feb 2010.
	Added -p n  switch to allow position output to be higher precision.  10 Nov 2012.
     Changed spot from int to long in the function re().  29 July 2013.  ( Thanks to Yuri D'Elia for this suggestion.)
      Changed function definitions, and misc other things to comply with c99
     requirements.  4 Mar 2014.
***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include "ms_new.h"

#define SITESINC 10

unsigned maxsites = SITESINC ;


//double *posit ;
//double segfac ;
//int count;
struct params pars ;
#pragma omp threadprivate(pars)

//#pragma omp threadprivate(segfac, count, pars)


int main_ms_ABLE(int ms_argc, char *ms_argv[], double **onetreePoisTable)
{
	int howmany;
	void getpars( int ms_argc, char *ms_argv[], int *howmany )  ;
	int gensam_ABLE(double **onetreePoisTable, int *crash_flag) ;
 	void freed2matrix(double **m, int x);
	void free_eventlist( struct devent *pt, int npop );

//	count=0;

	getpars( ms_argc, ms_argv, &howmany) ;

	int ms_thread_crash_flag = 0;
    gensam_ABLE(onetreePoisTable, &ms_thread_crash_flag);

	// based on Valgrind Memcheck
	free(pars.cp.config);
	free(pars.cp.size);
	free(pars.cp.alphag);
	free_eventlist(pars.cp.deventlist, pars.cp.npop);
	freed2matrix(pars.cp.mig_mat, pars.cp.npop);

	return ms_thread_crash_flag;
}


int gensam_ABLE(double **onetreePoisTable, int *crash_flag)
{
	int nsegs, i, j, k, seg, ns, start, end, len;
	struct segl *seglst, *segtre_mig(struct c_params *p, int *nsegs, int *crash_flag) ; /* used to be: [MAXSEG];  */
	int segsitesin,nsites;
	double theta;
	int nsam;
	void evalTreeBranchConfigs(struct node *ptree, int nsam, double *totbrlen);
	double ** d2matrix(int x, int y);
	void freed2matrix(double **m, int x);
	void storePoissonProbs(double **onetreeTable);


	nsites = pars.cp.nsites;
	seglst = segtre_mig(&(pars.cp), &nsegs, crash_flag);

	if (*crash_flag) {
		for (seg = 0, k = 0; k < nsegs; seg = seglst[seg].next, k++)
			free(seglst[seg].ptree);
		free(seglst);
		return 0;
	}

	nsam = pars.cp.nsam;
	segsitesin = pars.mp.segsitesin;
	theta = pars.mp.theta;

	ns = 0;

	int *onetreesegs = (int *) malloc((nsegs) * sizeof(int));
	double **totSegBrLen = d2matrix(nsegs, allBrClasses);
	double *totBrLen = (double *) malloc(brClass * sizeof(double));

	for (seg = 0, k = 0; k < nsegs; seg = seglst[seg].next, k++) {
		if ((pars.cp.r > 0.0) || (pars.cp.f > 0.0)) {
			end = (k < nsegs - 1 ? seglst[seglst[seg].next].beg - 1 : nsites - 1);
			start = seglst[seg].beg;
			len = end - start + 1;
			onetreesegs[k] = len;
		}

		evalTreeBranchConfigs(seglst[seg].ptree, nsam, totSegBrLen[k]);

		if ((segsitesin == 0) && (theta == 0.0) && (pars.mp.timeflag == 0))
			free(seglst[seg].ptree);
	}

	for (i = 1; i <= brClass; i++) {
		totBrLen[i - 1] = 0.0;

		for (k = 0; k < nsegs; k++) {
			double totFoldedSegBrLen = totSegBrLen[k][i - 1] + (foldBrClass * totSegBrLen[k][allBrClasses - i]);
			if (foldBrClass && ((i - 1) == (allBrClasses - i)))
				totFoldedSegBrLen /= 2;

			if ((pars.cp.r > 0.0) || (pars.cp.f > 0.0))
				totBrLen[i - 1] += totFoldedSegBrLen * onetreesegs[k];
			else
				totBrLen[i - 1] += totFoldedSegBrLen;
		}
		if ((pars.cp.r > 0.0) || (pars.cp.f > 0.0))
			totBrLen[i - 1] /= nsites;

//		printf("***Folded %d-ton branches***\n", i);
		if (totBrLen[i - 1] > 0.0) {
//			index j = mutClass-1 reserved for the marginal probabilities (i.e. gsl_cdf_poisson_Q())
			for (j = 0; j < mutClass - 1; j++) {
//				printf("%d : %5.5lf\n", j, gsl_ran_poisson_pdf(j,totBrLen[i-1]*pars.mp.theta));
				onetreePoisTable[i - 1][j] = gsl_ran_poisson_pdf(j, totBrLen[i - 1] * pars.mp.theta);
			}
			onetreePoisTable[i - 1][j] = gsl_cdf_poisson_Q(j - 1, totBrLen[i - 1] * pars.mp.theta);

//			printf(">%d : %5.5lf\n", j, gsl_cdf_poisson_Q(j,totBrLen[i-1]*pars.mp.theta));
//			printf("Total folded branch length = %5.5lf\n\n",totBrLen[i-1]);
		} else {
//			printf("Total folded branch length = 0\n\n");
			onetreePoisTable[i - 1][0] = 1.0;
			for (j = 1; j < mutClass; j++)
				onetreePoisTable[i - 1][j] = 0.0;
		}
	}

//	printf("\n");

	freed2matrix(totSegBrLen, nsegs);
	free(totBrLen);
	free(onetreesegs);

	// based on Valgrind Memcheck
	for (seg = 0, k = 0; k < nsegs; seg = seglst[seg].next, k++)
		free(seglst[seg].ptree);
	free(seglst);

	return (ns);
}





void ndes_and_branch_setup(struct node *ptree, int nsam) {
	int i ;

	for(i = 0; i < nsam; i++)
		ptree[i].ndes = 1 ;

	for(i = nsam; i < 2*nsam - 1; i++)
		ptree[i].ndes = 0 ;

	for(i = 0; i < 2*nsam - 2; i++) {
		ptree[(ptree[i].abv)].ndes += ptree[i].ndes ;
		ptree[i].brLength = ptree[ptree[i].abv].time - ptree[i].time;
	}
	//	ROOT!
	ptree[2*nsam -2].brLength = 0 ;
}

void evalTreeBranchConfigs(struct node *ptree, int nsam, double *totSegBrLenK) {
	int branch, pop, tip, ind, class;

	void ndes_and_branch_setup(struct node *ptree, int nsam);
	int ** d2int_matrix(int x, int y);
	void freed2int_matrix(int **m, int x);
	int getPopSampleStatus(int pop);
	int getBrConfigNum(int *brConfVec);

	int **branchConfig = d2int_matrix(2*nsam - 2, pars.cp.npop);
	int **sampledBranchConfig = d2int_matrix(2*nsam - 2, sampledPopsSize);

	ndes_and_branch_setup(ptree, nsam);

	for(branch = 0; branch < 2*nsam - 2; branch++)
		for(pop = 0; pop < pars.cp.npop; pop++)
			branchConfig[branch][pop] = 0;

	tip = 0;
	for(pop = 0; pop < pars.cp.npop; pop++) {
		for(ind = 1; ind <= pars.cp.config[pop]; ind++) {
			branchConfig[tip][pop] = 1;
			++tip;
		}
	}

	for(branch = 0; branch < 2*nsam - 2; branch++)
		for(pop = 0; pop < pars.cp.npop; pop++)
			if (ptree[branch].abv != 2*nsam - 2)
				branchConfig[ptree[branch].abv][pop] += branchConfig[branch][pop];

	for(class = 0; class < allBrClasses; class++)
		totSegBrLenK[class] = 0.0;

	//	if ghost/unsampled pops. have been specified
	if (sampledPopsSize != pars.cp.npop) {
		for(branch = 0; branch < 2*nsam - 2; branch++)
			for(pop = 0; pop < pars.cp.npop; pop++)
				if (getPopSampleStatus(pop))
					sampledBranchConfig[branch][pop] = branchConfig[branch][pop];

		for(branch = 0; branch < 2*nsam-2; branch++)
			totSegBrLenK[getBrConfigNum(sampledBranchConfig[branch])] += ptree[branch].brLength;
	}
	else
		for(branch = 0; branch < 2*nsam-2; branch++)
			totSegBrLenK[getBrConfigNum(branchConfig[branch])] += ptree[branch].brLength;

	//******************************************************************************
/*
	printf("\n");
	for(class = 0; class < allBrClasses; class++)
		printf("%5.3f\n",totSegBrLenK[class]);
	exit(-1);
*/
	//******************************************************************************

	freed2int_matrix(branchConfig, 2*nsam-2);
	freed2int_matrix(sampledBranchConfig, 2*nsam-2);
}




/*
int main_ms(int ms_argc, char *ms_argv[])
{

	int i, k, howmany, segsites, listX, listY;
	char **list, **cmatrix(), **tbsparamstrs ;
	FILE *pf, *fopen() ;
	char **list, **cmatrix(), **tbsparamstrs ;

	int i, howmany, listX;
	char **list, **cmatrix();
	FILE *fopen() ;
	double probss, tmrca, ttot ;
	void seedit( const char * ) ;
	void getpars( int ms_argc, char *ms_argv[], int *howmany )  ;
	int gensam( char **list, double *probss, double *ptmrca, double *pttot) ;
	void freecmatrix(char **m, int nsam);
 	void freed2matrix(double **m, int x);



	ntbs = 0 ;    these next few lines are for reading in parameters from a file (for each sample)
	tbsparamstrs = (char **)malloc( ms_argc*sizeof(char *) ) ;

	for( i=1; i<ms_argc; i++)
		printf("%s ",ms_argv[i]);
	printf("\n");
	for( i =0; i<ms_argc; i++) tbsparamstrs[i] = (char *)malloc(30*sizeof(char) ) ;
	for( i = 1; i<ms_argc ; i++)
			if( strcmp( ms_argv[i],"tbs") == 0 )  ms_argv[i] = tbsparamstrs[ ntbs++] ;


	count=0;

//	if( ntbs > 0 )  for( k=0; k<ntbs; k++)  scanf(" %s", tbsparamstrs[k] );
	getpars( ms_argc, ms_argv, &howmany) ;    results are stored in global variable, pars
//	if( !pars.commandlineseedflag ) seedit( "s");
//	pf = stdout ;

	if( pars.mp.segsitesin ==  0 ) {
	     list = cmatrix(pars.cp.nsam,maxsites+1);
//	     posit = (double *)malloc( (unsigned)( maxsites*sizeof( double)) ) ;
	     listX = pars.cp.nsam;
//	     listY = maxsites+1;
	}
	else {
	     list = cmatrix(pars.cp.nsam, pars.mp.segsitesin+1 ) ;
//	     posit = (double *)malloc( (unsigned)( pars.mp.segsitesin*sizeof( double)) ) ;
	     if( pars.mp.theta > 0.0 ){
		    segfac = 1.0 ;
		    for(  i= pars.mp.segsitesin; i > 1; i--) segfac *= i ;
		 }
	     listX = pars.cp.nsam;
//	     listY = pars.mp.segsitesin+1;
	}

    while( howmany-count++ ) {

	   if( (ntbs > 0) && (count >1 ) ){
	         for( k=0; k<ntbs; k++){
			    if( scanf(" %s", tbsparamstrs[k]) == EOF ){
			       if( !pars.commandlineseedflag ) seedit( "end" );
				   exit(0);
				}
			 }
			 getpars( ms_argc, ms_argv, &howmany) ;
	   }

       fprintf(pf,"\n//");
	   if( ntbs >0 ){
			for(k=0; k< ntbs; k++) printf("\t%s", tbsparamstrs[k] ) ;
		}
		printf("\n");

        segsites = gensam( list, &probss, &tmrca, &ttot) ;

        gensam( list, &probss, &tmrca, &ttot);

        if (ms_crash_flag)
        	break;


  		if( pars.mp.timeflag ) fprintf(pf,"time:\t%lf\t%lf\n",tmrca, ttot ) ;
        if( (segsites > 0 ) || ( pars.mp.theta > 0.0 ) ) {
   	       if( (pars.mp.segsitesin > 0 ) && ( pars.mp.theta > 0.0 ))
		       fprintf(pf,"prob: %g\n", probss ) ;
           fprintf(pf,"segsites: %d\n",segsites);
		   if( segsites > 0 )	fprintf(pf,"positions: ");
		   for( i=0; i<segsites; i++)
              fprintf(pf,"%6.*lf ", pars.output_precision,posit[i] );
           fprintf(pf,"\n");
	       if( segsites > 0 )
		          for(i=0;i<pars.cp.nsam; i++) { fprintf(pf,"%s\n", list[i] ); }
	    }

    }
//	if( !pars.commandlineseedflag ) seedit( "end" );

	freecmatrix(list, listX);

	// based on Valgrind Memcheck
	free(pars.cp.config);
	free(pars.cp.size);
	free(pars.cp.alphag);
	free(pars.cp.deventlist);
	freed2matrix(pars.cp.mig_mat, pars.cp.npop);

	return 0;
}



	int 
gensam( char **list, double *pprobss, double *ptmrca, double *pttot)
{

	int nsegs, i, j, k, seg, ns, start, end, len, segsit ;
	struct segl *seglst, *segtre_mig(struct c_params *p, int *nsegs ) ;  used to be: [MAXSEG];
	double nsinv,  tseg, tt, ttime(struct node *, int nsam), ttimemf(struct node *, int nsam, int mfreq) ;
	double *pk;
	int *ss;
	int segsitesin,nsites;
	double theta, es ;
	int nsam, mfreq ;

	int nsegs, i, j, k, seg, ns, start, end, len;
	struct segl *seglst, *segtre_mig(struct c_params *p, int *nsegs ) ;  used to be: [MAXSEG];
	double ttime(struct node *, int nsam), ttimemf(struct node *, int nsam, int mfreq) ;
	int segsitesin,nsites;
	double theta;
	int nsam;
	void prtree( struct node *ptree, int nsam);
	void make_gametes(int nsam, int mfreq,  struct node *ptree, double tt, int newsites, int ns, char **list );
 	void ndes_setup( struct node *, int nsam );
 	void evalTreeBranchConfigs(struct node *ptree, int nsam, double *totbrlen);
 	double ** d2matrix(int x, int y);
 	void freed2matrix(double **m, int x);
 	void storePoissonProbs(double **onetreeTable);


	nsites = pars.cp.nsites ;
//	nsinv = 1./nsites;
	seglst = segtre_mig(&(pars.cp),  &nsegs ) ;
	
	if (ms_crash_flag) {
	    for(seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++)
	    	free(seglst[seg].ptree);
		return 0;
	}

	nsam = pars.cp.nsam;
	segsitesin = pars.mp.segsitesin ;
	theta = pars.mp.theta ;
//	mfreq = pars.mp.mfreq ;

//	if( pars.mp.treeflag ) {
	  	ns = 0 ;


	  	double **onetreeTable = d2matrix(brClass, mutClass);
	  	int *onetreesegs = (int *) malloc((nsegs) * sizeof(int));
	  	double **totSegBrLen = d2matrix(nsegs, allBrClasses);
	  	double *totBrLen = (double *) malloc(brClass * sizeof(double));

	    for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) {
	    	if( (pars.cp.r > 0.0 ) || (pars.cp.f > 0.0) ){
	    		end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
	    		start = seglst[seg].beg ;
	    		len = end - start + 1 ;
//	    		fprintf(stdout,"[%d]", len);
	    		onetreesegs[k] = len;
	    	}
//	    	prtree( seglst[seg].ptree, nsam) ;

	    	evalTreeBranchConfigs(seglst[seg].ptree, nsam, totSegBrLen[k]);

	    	if( (segsitesin == 0) && ( theta == 0.0 ) && ( pars.mp.timeflag == 0 ) )
	    		free(seglst[seg].ptree) ;
	    }

	    for(i = 1; i <= brClass; i++) {
	    	totBrLen[i-1] = 0.0;

		    for(k = 0; k < nsegs; k++) {
		    	double totFoldedSegBrLen = totSegBrLen[k][i-1] + (foldBrClass * totSegBrLen[k][allBrClasses-i]);
			    if (foldBrClass && ((i-1) == (allBrClasses-i)))
			    	totFoldedSegBrLen /= 2;

		    	if( (pars.cp.r > 0.0 ) || (pars.cp.f > 0.0) )
		    		totBrLen[i-1] += totFoldedSegBrLen * onetreesegs[k];
		    	else
		    		totBrLen[i-1] += totFoldedSegBrLen;
		    }
	    	if( (pars.cp.r > 0.0 ) || (pars.cp.f > 0.0) )
	    		totBrLen[i-1] /= nsites;

//	    	printf("***Folded %d-ton branches***\n", i);
    		if (totBrLen[i-1] > 0.0) {
    			// index j = mutClass-1 reserved for the marginal probabilities (i.e. gsl_cdf_poisson_Q())
    			for(j = 0; j < mutClass-1; j++) {
//    				printf("%d : %5.5lf\n", j, gsl_ran_poisson_pdf(j,totBrLen[i-1]*pars.mp.theta));
    				onetreeTable[i-1][j] = gsl_ran_poisson_pdf(j,totBrLen[i-1]*pars.mp.theta);
    			}
    			onetreeTable[i-1][j] = gsl_cdf_poisson_Q(j-1,totBrLen[i-1]*pars.mp.theta);

//    			printf(">%d : %5.5lf\n", j, gsl_cdf_poisson_Q(j,totBrLen[i-1]*pars.mp.theta));
//    			printf("Total folded branch length = %5.5lf\n\n",totBrLen[i-1]);
    		}
    		else {
//    			printf("Total folded branch length = 0\n\n");
    			onetreeTable[i-1][0] = 1.0;
    			for(j = 1; j < mutClass; j++)
    				onetreeTable[i-1][j] = 0.0;
    		}
    	}

	    storePoissonProbs(onetreeTable);

//	    printf("\n");


	    freed2matrix(totSegBrLen, nsegs);
	    free(totBrLen);
	    free(onetreesegs);

		// based on Valgrind Memcheck
	    for(seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++)
	    	free(seglst[seg].ptree);

//	}


	if( pars.mp.timeflag ) {
      tt = 0.0 ;
	  for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) {
		if( mfreq > 1 ) ndes_setup( seglst[seg].ptree, nsam );
		end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
		start = seglst[seg].beg ;
		if( (nsegs==1) || ( ( start <= nsites/2) && ( end >= nsites/2 ) ) )
		  *ptmrca = (seglst[seg].ptree + 2*nsam-2) -> time ;
		len = end - start + 1 ;
		tseg = len/(double)nsites ;
		if( mfreq == 1 ) tt += ttime(seglst[seg].ptree,nsam)*tseg ;
		else tt += ttimemf(seglst[seg].ptree,nsam, mfreq)*tseg ;
		if( (segsitesin == 0) && ( theta == 0.0 )  )
	  	      free(seglst[seg].ptree) ;
	    }
		*pttot = tt ;
	 }

    if( (segsitesin == 0) && ( theta > 0.0)   ) {
	  ns = 0 ;
	  for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) {
		if( mfreq > 1 ) ndes_setup( seglst[seg].ptree, nsam );
		end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
		start = seglst[seg].beg ;
		len = end - start + 1 ;
		tseg = len*(theta/nsites) ;
		if( mfreq == 1) tt = ttime(seglst[seg].ptree, nsam);
                else tt = ttimemf(seglst[seg].ptree, nsam, mfreq );
		segsit = poisso( tseg*tt );
		if( (segsit + ns) >= maxsites ) {
			maxsites = segsit + ns + SITESINC ;
			posit = (double *)realloc(posit, maxsites*sizeof(double) ) ;
			  biggerlist(nsam, list) ;
		}
		make_gametes(nsam,mfreq,seglst[seg].ptree,tt, segsit, ns, list );
		free(seglst[seg].ptree) ;
		locate(segsit,start*nsinv, len*nsinv,posit+ns);
		ns += segsit;
	  }
    }
   else if( segsitesin > 0 ) {

        pk = (double *)malloc((unsigned)(nsegs*sizeof(double)));
        ss = (int *)malloc((unsigned)(nsegs*sizeof(int)));
        if( (pk==NULL) || (ss==NULL) ) perror("malloc error. gensam.2");


	  tt = 0.0 ;
	  for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) {
		if( mfreq > 1 ) ndes_setup( seglst[seg].ptree, nsam );
		end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
		start = seglst[seg].beg ;
		len = end - start + 1 ;
		tseg = len/(double)nsites ;
               if( mfreq == 1 ) pk[k] = ttime(seglst[seg].ptree,nsam)*tseg ;
               else pk[k] = ttimemf(seglst[seg].ptree,nsam, mfreq)*tseg ;
                 tt += pk[k] ;
	  }
	  if( theta > 0.0 ) {
	    es = theta * tt ;
	    *pprobss = exp( -es )*pow( es, (double) segsitesin) / segfac ;
	  }
	  if( tt > 0.0 ) {
          for (k=0;k<nsegs;k++) pk[k] /= tt ;
          mnmial(segsitesin,nsegs,pk,ss);
	  }
	  else
            for( k=0; k<nsegs; k++) ss[k] = 0 ;
	  ns = 0 ;
	  for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) {
		 end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
		 start = seglst[seg].beg ;
		 len = end - start + 1 ;
		 tseg = len/(double)nsites;
		 make_gametes(nsam,mfreq,seglst[seg].ptree,tt*pk[k]/tseg, ss[k], ns, list);

		 free(seglst[seg].ptree) ;
		 locate(ss[k],start*nsinv, len*nsinv,posit+ns);
		 ns += ss[k] ;
	  }
	  free(pk);
	  free(ss);

    }
	for(i=0;i<nsam;i++) list[i][ns] = '\0' ;

	return( ns ) ;
}

void ndes_setup(struct node *ptree, int nsam )
{
	int i ;

	for( i=0; i<nsam; i++) (ptree+i)->ndes = 1 ;
	for( i = nsam; i< 2*nsam -1; i++) (ptree+i)->ndes = 0 ;
	for( i= 0; i< 2*nsam -2 ; i++)  (ptree+((ptree+i)->abv))->ndes += (ptree+i)->ndes ;

}
*/




/*
	void
biggerlist(int nsam,  char **list )
{
	int i;

  fprintf(stderr,"maxsites: %d\n",maxsites);
	for( i=0; i<nsam; i++){
	   list[i] = (char *)realloc( list[i],maxsites*sizeof(char) ) ;
	   if( list[i] == NULL ) perror( "realloc error. bigger");
	   }
}
	   


 allocates space for gametes (character strings)
	char **
cmatrix(nsam,len)
	int nsam, len;
{
	int i;
	char **m;

	if( ! ( m = (char **) malloc( (unsigned) nsam*sizeof( char* ) ) ) )
	   perror("alloc error in cmatrix") ;
	for( i=0; i<nsam; i++) {
		if( ! ( m[i] = (char *) malloc( (unsigned) len*sizeof( char ) )))
			perror("alloc error in cmatrix 2");
		}
	return( m );
}



void freecmatrix(char **m, int nsam)
{
	int i;
	for( i=nsam-1; i>=0; i--)
		free(m[i]);
	free(m);
}
*/



/* allocates space for mutation probabilities on branch classes (3-dim double matrix) */
double *** d3matrix(int x, int y, int z)
{
	int i, j;
	double ***tab;

	if( ! ( tab = (double ***) malloc( (unsigned) x*sizeof( double** ) ) ) )
		perror("alloc error in d3matrix") ;
	for( i=0; i<x; i++) {
		if( ! ( tab[i] = (double **) malloc( (unsigned) y*sizeof( double* ) )))
			perror("alloc error in d3matrix. 2");
		for( j=0; j<y; j++) {
			if( ! ( tab[i][j] = (double *) malloc( (unsigned) z*sizeof( double ) )))
				perror("alloc error in d3matrix 3");
		}
	}

	return( tab );
}


void freed3matrix(double ***m, int x, int y)
{
	int i, j;
	for( i=x-1; i>=0; i--) {
		for( j=y-1; j>=0; j--)
			free(m[i][j]);
		free(m[i]);
	}
	free(m);
}




/* allocates space for mutation probabilities on branch classes (2-dim double matrix) */
double ** d2matrix(int x, int y)
{
	int i;
	double **tab;

	if( ! ( tab = (double **) malloc( (unsigned) x*sizeof( double* ) ) ) )
		perror("alloc error in d2matrix") ;
	for( i=0; i<x; i++) {
		if( ! ( tab[i] = (double *) malloc( (unsigned) y*sizeof( double ) )))
			perror("alloc error in d2matrix. 2");
	}

	return( tab );
}


void freed2matrix(double **m, int x)
{
	int i;
	for( i=x-1; i>=0; i--)
		free(m[i]);
	free(m);
}




/* allocates sfreed2matrixpace for mutation probabilities on branch classes (2-dim int matrix) */
int ** d2int_matrix(int x, int y)
{
	int i;
	int **tab;

	if( ! ( tab = (int **) malloc( (unsigned) x*sizeof( int* ) ) ) )
		perror("alloc error in dmatrix") ;
	for( i=0; i<x; i++) {
		if( ! ( tab[i] = (int *) malloc( (unsigned) y*sizeof( int ) )))
			perror("alloc error in dmatrix. 2");
	}

	return( tab );
}


void freed2int_matrix(int **m, int x)
{
	int i;
	for( i=x-1; i>=0; i--)
		free(m[i]);
	free(m);
}




/*
	void
locate(int n,double beg, double len,double *ptr)
{
	int i;

	ordran(n,ptr);
	for(i=0; i<n; i++)
		ptr[i] = beg + ptr[i]*len ;

}
*/

int NSEEDS = 3 ;

  void
getpars(int ms_argc, char *ms_argv[], int *phowmany )
{
//	int arg, i, j, sum , pop , argstart, npop , npop2, pop2 ;
	int arg, i, j, sum , pop , npop , npop2, pop2 ;
	double migr, mij, psize, palpha ;
	void addtoelist( struct devent *pt, struct devent *elist ); 
	void argcheck( int arg, int ms_argc, char ** ) ;
	int commandlineseed( char ** ) ;
	void free_eventlist( struct devent *pt, int npop );
	struct devent *ptemp , *pt ;
//	FILE *pf ;
	char ch3 ;
	

//  if( count == 0 ) {
	if( ms_argc < 4 ){ fprintf(stderr,"Too few command line arguments\n"); usage();}
	pars.cp.nsam = atoi( ms_argv[1] );
	if( pars.cp.nsam <= 0 ) { fprintf(stderr,"First argument error. nsam <= 0. \n"); usage();}
	*phowmany = atoi( ms_argv[2] );
	if( *phowmany  <= 0 ) { fprintf(stderr,"Second argument error. howmany <= 0. \n"); usage();}
//	pars.commandlineseedflag = 0 ;
	  pars.output_precision = 4 ;
	pars.cp.r = pars.mp.theta =  pars.cp.f = 0.0 ;
	pars.cp.track_len = 0. ;
	pars.cp.npop = npop = 1 ;
	pars.cp.mig_mat = (double **)malloc( (unsigned) sizeof( double *) );
	pars.cp.mig_mat[0] = (double *)malloc( (unsigned)sizeof(double ));
	pars.cp.mig_mat[0][0] =  0.0 ;
	pars.mp.segsitesin = 0 ;
	pars.mp.treeflag = 0 ;
 	pars.mp.timeflag = 0 ;
       pars.mp.mfreq = 1 ;
	pars.cp.config = (int *) malloc( (unsigned)(( pars.cp.npop +1 ) *sizeof( int)) );
	(pars.cp.config)[0] = pars.cp.nsam ;
	pars.cp.size= (double *) malloc( (unsigned)( pars.cp.npop *sizeof( double )) );
	(pars.cp.size)[0] = 1.0  ;
	pars.cp.alphag = (double *) malloc( (unsigned)(( pars.cp.npop ) *sizeof( double )) );
	(pars.cp.alphag)[0] = 0.0  ;
	pars.cp.nsites = 2 ;
/*
  }
  else{
	npop = pars.cp.npop ;
	free_eventlist( pars.cp.deventlist, npop );
  }
*/
  	pars.cp.deventlist = NULL ;

	arg = 3 ;

	while( arg < (ms_argc - 2) ){
		if( ms_argv[arg][0] != '-' ) { fprintf(stderr," argument should be -%s ?\n", ms_argv[arg]); usage();}
		switch ( ms_argv[arg][1] ){
/*
			case 'f' :
//				if( ntbs > 0 ) { fprintf(stderr," can't use tbs args and -f option.\n"); exit(1); }
				arg++;
				argcheck( arg, ms_argc, ms_argv);
				pf = fopen( ms_argv[arg], "r" ) ;
				if( pf == NULL ) {fprintf(stderr," no parameter file %s\n", ms_argv[arg] ); exit(0);}
				arg++;
				ms_argc++ ;
				ms_argv = (char **)malloc(  (unsigned)(ms_argc+1)*sizeof( char *) ) ;
				ms_argv[arg] = (char *)malloc( (unsigned)(20*sizeof( char )) ) ;
				argstart = arg ;
				while( fscanf(pf," %s", ms_argv[arg]) != EOF ) {
					arg++;
					ms_argc++;
					ms_argv = (char **)realloc( ms_argv, (unsigned)ms_argc*sizeof( char*) ) ;
				        ms_argv[arg] = (char *)malloc( (unsigned)(20*sizeof( char )) ) ;
					}
				fclose(pf);
				ms_argc--;
				arg = argstart ;
				break;
*/
			case 'r' : 
				arg++;
				argcheck( arg, ms_argc, ms_argv);
				pars.cp.r = atof(  ms_argv[arg++] );
				argcheck( arg, ms_argc, ms_argv);
				pars.cp.nsites = atoi( ms_argv[arg++]);
				if( pars.cp.nsites <2 ){
					fprintf(stderr,"with -r option must specify both rec_rate and nsites>1\n");
					usage();
					}
				break;	
/*
			case 'p' :
				arg++;
				argcheck(arg,ms_argc,ms_argv);
				pars.output_precision = atoi( ms_argv[arg++] ) ;
				break;
			case 'c' : 
				arg++;
				argcheck( arg, ms_argc, ms_argv);
				pars.cp.f = atof(  ms_argv[arg++] );
				argcheck( arg, ms_argc, ms_argv);
				pars.cp.track_len = atof( ms_argv[arg++]);
				if( pars.cp.track_len <1. ){
					fprintf(stderr,"with -c option must specify both f and track_len>0\n");
					usage();
					}
				break;		
*/
			case 't' : 
				arg++;
				argcheck( arg, ms_argc, ms_argv);
				pars.mp.theta = atof(  ms_argv[arg++] );
				break;
/*
			case 's' : 
				arg++;
				argcheck( arg, ms_argc, ms_argv);
				if( ms_argv[arg-1][2] == 'e' ){   command line seeds
//					pars.commandlineseedflag = 1 ;
//					if( count == 0 ) nseeds = commandlineseed(ms_argv+arg );
					arg += nseeds ;
				}
				else {
				    pars.mp.segsitesin = atoi(  ms_argv[arg++] );
				}
				break;
			case 'F' : 
				arg++;
				argcheck( arg, ms_argc, ms_argv);
				pars.mp.mfreq = atoi(  ms_argv[arg++] );
                                if( (pars.mp.mfreq < 2 ) || (pars.mp.mfreq > pars.cp.nsam/2 ) ){
                                    fprintf(stderr," mfreq must be >= 2 and <= nsam/2.\n");
                                    usage();
                                    }
				break;
*/
			case 'T' : 
				pars.mp.treeflag = 1 ;
				arg++;
				break;
/*
			case 'L' : 
				pars.mp.timeflag = 1 ;
				arg++;
				break;
*/
			case 'I' : 
			    arg++;
//			    if( count == 0 ) {
				argcheck( arg, ms_argc, ms_argv);
			       	pars.cp.npop = atoi( ms_argv[arg]);
			        pars.cp.config = (int *) realloc( pars.cp.config, (unsigned)( pars.cp.npop*sizeof( int)));
				npop = pars.cp.npop ;
//				}
			    arg++;
			    for( i=0; i< pars.cp.npop; i++) {
				argcheck( arg, ms_argc, ms_argv);
				pars.cp.config[i] = atoi( ms_argv[arg++]);
				}
//			    if( count == 0 ){
				pars.cp.mig_mat = 
                                        (double **)realloc(pars.cp.mig_mat, (unsigned)(pars.cp.npop*sizeof(double *) )) ;
				pars.cp.mig_mat[0] = 
                                         (double *)realloc(pars.cp.mig_mat[0], (unsigned)( pars.cp.npop*sizeof(double)));
				for(i=1; i<pars.cp.npop; i++) pars.cp.mig_mat[i] = 
                                         (double *)malloc( (unsigned)( pars.cp.npop*sizeof(double)));
				pars.cp.size = (double *)realloc( pars.cp.size, (unsigned)( pars.cp.npop*sizeof( double )));
				pars.cp.alphag = 
                                          (double *) realloc( pars.cp.alphag, (unsigned)( pars.cp.npop*sizeof( double )));
			        for( i=1; i< pars.cp.npop ; i++) {
				   (pars.cp.size)[i] = (pars.cp.size)[0]  ;
				   (pars.cp.alphag)[i] = (pars.cp.alphag)[0] ;
				   }
//			        }
			     if( (arg <ms_argc) && ( ms_argv[arg][0] != '-' ) ) {
				argcheck( arg, ms_argc, ms_argv);
				migr = atof(  ms_argv[arg++] );
				}
			     else migr = 0.0 ;
			     for( i=0; i<pars.cp.npop; i++) 
				    for( j=0; j<pars.cp.npop; j++) pars.cp.mig_mat[i][j] = migr/(pars.cp.npop-1) ;
			     for( i=0; i< pars.cp.npop; i++) pars.cp.mig_mat[i][i] = migr ;
			     break;
			case 'm' :
			     if( npop < 2 ) { fprintf(stderr,"Must use -I option first.\n"); usage();}
			     if( ms_argv[arg][2] == 'a' ) {
				    arg++;
				    for( pop = 0; pop <npop; pop++)
				      for( pop2 = 0; pop2 <npop; pop2++){
					     argcheck( arg, ms_argc, ms_argv);
					     pars.cp.mig_mat[pop][pop2]= atof( ms_argv[arg++] ) ;
					  }
				    for( pop = 0; pop < npop; pop++) {
					  pars.cp.mig_mat[pop][pop] = 0.0 ;
					  for( pop2 = 0; pop2 < npop; pop2++){
					    if( pop2 != pop ) pars.cp.mig_mat[pop][pop] += pars.cp.mig_mat[pop][pop2] ;
					  }
				    }	
				}
			    else {
		             arg++;
			         argcheck( arg, ms_argc, ms_argv);
		             i = atoi( ms_argv[arg++] ) -1;
			         argcheck( arg, ms_argc, ms_argv);
		             j = atoi( ms_argv[arg++] ) -1;
			         argcheck( arg, ms_argc, ms_argv);
		             mij = atof( ms_argv[arg++] );
		             pars.cp.mig_mat[i][i] += mij -  pars.cp.mig_mat[i][j]  ;
		             pars.cp.mig_mat[i][j] = mij;
			    }
				break;
			case 'n' :
			     if( npop < 2 ) { fprintf(stderr,"Must use -I option first.\n"); usage();}
			    arg++;
			    argcheck( arg, ms_argc, ms_argv);
			    pop = atoi( ms_argv[arg++] ) -1;
			    argcheck( arg, ms_argc, ms_argv);
			    psize = atof( ms_argv[arg++] );
			    pars.cp.size[pop] = psize ;
			   break;
			case 'g' :
			     if( npop < 2 ) { fprintf(stderr,"Must use -I option first.\n"); usage();}
			    arg++;
			    argcheck( arg, ms_argc, ms_argv);
			    pop = atoi( ms_argv[arg++] ) -1;
			    if( arg >= ms_argc ) { fprintf(stderr,"Not enough arg's after -G.\n"); usage(); }
			    palpha = atof( ms_argv[arg++] );
			    pars.cp.alphag[pop] = palpha ;
			   break;
			case 'G' :
			    arg++;
			    if( arg >= ms_argc ) { fprintf(stderr,"Not enough arg's after -G.\n"); usage(); }
			    palpha = atof( ms_argv[arg++] );
			    for( i=0; i<pars.cp.npop; i++) 
			       pars.cp.alphag[i] = palpha ;
			   break;
			case 'e' :
			    pt = (struct devent *)malloc( sizeof( struct devent) ) ;
			    pt->detype = ms_argv[arg][2] ;
			    ch3 = ms_argv[arg][3] ;
			    arg++;
			    argcheck( arg, ms_argc, ms_argv);
			    pt->time = atof( ms_argv[arg++] ) ;
			    pt->nextde = NULL ;
			    if( pars.cp.deventlist == NULL ) 
				    pars.cp.deventlist = pt ;
			    else if ( pt->time < pars.cp.deventlist->time ) { 
				    ptemp = pars.cp.deventlist ;
				    pars.cp.deventlist = pt ;
				    pt->nextde = ptemp ;	
				}	
			    else
				   addtoelist( pt, pars.cp.deventlist ) ;
			    switch( pt->detype ) {
				case 'N' :
			          argcheck( arg, ms_argc, ms_argv);
				      pt->paramv = atof( ms_argv[arg++] ) ;
				      break;
				case 'G' :
				  if( arg >= ms_argc ) { fprintf(stderr,"Not enough arg's after -eG.\n"); usage(); }
				  pt->paramv = atof( ms_argv[arg++] ) ;
				  break;
				case 'M' :
				    argcheck( arg, ms_argc, ms_argv);
				    pt->paramv = atof( ms_argv[arg++] ) ;
				    break;
				case 'n' :
			          argcheck( arg, ms_argc, ms_argv);
				  pt->popi = atoi( ms_argv[arg++] ) -1 ;
			          argcheck( arg, ms_argc, ms_argv);
				  pt->paramv = atof( ms_argv[arg++] ) ;
				  break;
				case 'g' :
			          argcheck( arg, ms_argc, ms_argv);
				  pt->popi = atoi( ms_argv[arg++] ) -1 ;
				  if( arg >= ms_argc ) { fprintf(stderr,"Not enough arg's after -eg.\n"); usage(); }
				  pt->paramv = atof( ms_argv[arg++] ) ;
				  break;
				case 's' :
			          argcheck( arg, ms_argc, ms_argv);
				  pt->popi = atoi( ms_argv[arg++] ) -1 ;
			          argcheck( arg, ms_argc, ms_argv);
				  pt->paramv = atof( ms_argv[arg++] ) ;
				  break;
				case 'm' :
				  if( ch3 == 'a' ) {
				     pt->detype = 'a' ;
				     argcheck( arg, ms_argc, ms_argv);
				     npop2 = atoi( ms_argv[arg++] ) ;
				     pt->mat = (double **)malloc( (unsigned)npop2*sizeof( double *) ) ;
				     for( pop =0; pop <npop2; pop++){
					   (pt->mat)[pop] = (double *)malloc( (unsigned)npop2*sizeof( double) );
					   for( i=0; i<npop2; i++){
					     if( i == pop ) arg++;
					     else {
				               argcheck( arg, ms_argc, ms_argv);
					       (pt->mat)[pop][i] = atof( ms_argv[arg++] ) ;
					     }
					   }
				     }
				     for( pop = 0; pop < npop2; pop++) {
					    (pt->mat)[pop][pop] = 0.0 ;
					    for( pop2 = 0; pop2 < npop2; pop2++){
					       if( pop2 != pop ) (pt->mat)[pop][pop] += (pt->mat)[pop][pop2] ;
					    }
				     }	
				  }
				  else {
			            argcheck( arg, ms_argc, ms_argv);
				        pt->popi = atoi( ms_argv[arg++] ) -1 ;
			            argcheck( arg, ms_argc, ms_argv);
				        pt->popj = atoi( ms_argv[arg++] ) -1 ;
			            argcheck( arg, ms_argc, ms_argv);
				        pt->paramv = atof( ms_argv[arg++] ) ;
				  }
				  break;
				case 'j' :
			          argcheck( arg, ms_argc, ms_argv);
				  pt->popi = atoi( ms_argv[arg++] ) -1 ;
			          argcheck( arg, ms_argc, ms_argv);
				  pt->popj = atoi( ms_argv[arg++] ) -1 ;
				  break;
				default: fprintf(stderr,"e event\n");  usage();
			    }
			 break;
			default: fprintf(stderr," option default\n");  usage() ;
			}
		}
		if( (pars.mp.theta == 0.0) && ( pars.mp.segsitesin == 0 ) && ( pars.mp.treeflag == 0 ) && (pars.mp.timeflag == 0) ) {
			fprintf(stderr," either -s or -t or -T option must be used. \n");
			usage();
			exit(1);
			}
		sum = 0 ;
		for( i=0; i< pars.cp.npop; i++) sum += (pars.cp.config)[i] ;
		if( sum != pars.cp.nsam ) {
			fprintf(stderr," sum sample sizes != nsam\n");
			usage();
			exit(1);
			}
}


	void
argcheck( int arg, int ms_argc, char *ms_argv[] )
{
	if( (arg >= ms_argc ) || ( ms_argv[arg][0] == '-') ) {
	   fprintf(stderr,"not enough arguments after %s\n", ms_argv[arg-1] ) ;
	   fprintf(stderr,"For usage type: ms<return>\n");
	   exit(0);
	  }
}
	
	void
usage()
{
fprintf(stderr,"usage: ms nsam howmany \n");
fprintf(stderr,"  Options: \n"); 
fprintf(stderr,"\t -t theta   (this option and/or the next must be used. Theta = 4*N0*u )\n");
//fprintf(stderr,"\t -s segsites   ( fixed number of segregating sites)\n");
fprintf(stderr,"\t -T          (Output gene tree.)\n");
//fprintf(stderr,"\t -F minfreq     Output only sites with freq of minor allele >= minfreq.\n");
fprintf(stderr,"\t -r rho nsites     (rho here is 4Nc)\n");
//fprintf(stderr,"\t\t -c f track_len   (f = ratio of conversion rate to rec rate. tracklen is mean length.) \n");
//fprintf(stderr,"\t\t\t if rho = 0.,  f = 4*N0*g, with g the gene conversion rate.\n");
fprintf(stderr,"\t -G alpha  ( N(t) = N0*exp(-alpha*t) .  alpha = -log(Np/Nr)/t\n");      
fprintf(stderr,"\t -I npop n1 n2 ... [mig_rate] (all elements of mig matrix set to mig_rate/(npop-1) \n");    
fprintf(stderr,"\t\t -m i j m_ij    (i,j-th element of mig matrix set to m_ij.)\n"); 
fprintf(stderr,"\t\t -ma m_11 m_12 m_13 m_21 m_22 m_23 ...(Assign values to elements of migration matrix.)\n"); 
fprintf(stderr,"\t\t -n i size_i   (popi has size set to size_i*N0 \n");
fprintf(stderr,"\t\t -g i alpha_i  (If used must appear after -M option.)\n"); 
fprintf(stderr,"\t   The following options modify parameters at the time 't' specified as the first argument:\n");
fprintf(stderr,"\t -eG t alpha  (Modify growth rate of all pop's.)\n");     
fprintf(stderr,"\t -eg t i alpha_i  (Modify growth rate of pop i.) \n");    
fprintf(stderr,"\t -eM t mig_rate   (Modify the mig matrix so all elements are mig_rate/(npop-1)\n"); 
fprintf(stderr,"\t -em t i j m_ij    (i,j-th element of mig matrix set to m_ij at time t )\n"); 
fprintf(stderr,"\t -ema t npop  m_11 m_12 m_13 m_21 m_22 m_23 ...(Assign values to elements of migration matrix.)\n");  
fprintf(stderr,"\t -eN t size  (Modify pop sizes. New sizes = size*N0 ) \n");    
fprintf(stderr,"\t -en t i size_i  (Modify pop size of pop i.  New size of popi = size_i*N0 .)\n");
fprintf(stderr,"\t -es t i proportion  (Split: pop i -> pop-i + pop-npop, npop increases by 1.\n");    
fprintf(stderr,"\t\t proportion is probability that each lineage stays in pop-i. (p, 1-p are admixt. proport.\n");
fprintf(stderr,"\t\t Size of pop npop is set to N0 and alpha = 0.0 , size and alpha of pop i are unchanged.\n");
fprintf(stderr,"\t -ej t i j   ( Join lineages in pop i and pop j into pop j\n");
fprintf(stderr,"\t\t  size, alpha and M are unchanged.\n");  
//fprintf(stderr,"\t  -f filename     ( Read command line arguments from file filename.)\n");
//fprintf(stderr,"\t  -p n ( Specifies the precision of the position output.  n is the number of digits after the decimal.)\n");
fprintf(stderr," See msdoc.pdf for explanation of these parameters.\n");

exit(1);
}
	void
addtoelist( struct devent *pt, struct devent *elist ) 
{
//	struct devent *plast, *pevent, *ptemp  ;
	struct devent *plast = NULL, *pevent, *ptemp  ;

	pevent = elist ;
	while(  (pevent != NULL ) && ( pevent->time <= pt->time ) )  {
		plast = pevent ;
		pevent = pevent->nextde ;
		}
	ptemp = plast->nextde ;
	plast->nextde = pt ;
	pt->nextde = ptemp ;
}

	void 
free_eventlist( struct devent *pt, int npop )
{
   struct devent *next ;
   int pop ;
   
   while( pt != NULL){
	  next = pt->nextde ;
	  if( pt->detype == 'a' ) {
	     for( pop = 0; pop < npop; pop++) free( (pt->mat)[pop] );
		 free( pt->mat );
	  }
	  free(pt);
	  pt = next ;
   }
}

	
/************ make_gametes.c  *******************************************
*
*
*****************************************************************************/

/*
#define STATE1 '1'
#define STATE2 '0'

	void
make_gametes(int nsam, int mfreq, struct node *ptree, double tt, int newsites, int ns, char **list )
{
	int  tip, j,  node ;
        int pickb(int nsam, struct node *ptree, double tt), 
            pickbmf(int nsam, int mfreq, struct node *ptree, double tt) ;

	for(  j=ns; j< ns+newsites ;  j++ ) {
		if( mfreq == 1 ) node = pickb(  nsam, ptree, tt);
		else node = pickbmf(  nsam, mfreq, ptree, tt);
		for( tip=0; tip < nsam ; tip++) {
		   if( tdesn(ptree, tip, node) ) list[tip][j] = STATE1 ;
		   else list[tip][j] = STATE2 ;
		   }
		}
}
*/


/***  ttime.c : Returns the total time in the tree, *ptree, with nsam tips. **/

/*
	double
ttime( ptree, nsam)
	struct node *ptree;
	int nsam;
{
	double t;
	int i;

	t = (ptree + 2*nsam-2) -> time ;
	for( i=nsam; i< 2*nsam-1 ; i++)
		t += (ptree + i)-> time ;
	return(t);
}


	double
ttimemf( ptree, nsam, mfreq)
	struct node *ptree;
	int nsam, mfreq;
{
	double t;
	int i;

	t = 0. ;
	for( i=0;  i< 2*nsam-2  ; i++)
	  if( ( (ptree+i)->ndes >= mfreq )  && ( (ptree+i)->ndes <= nsam-mfreq) )
		t += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
	return(t);
}


	void
prtree( ptree, nsam)
	struct node *ptree;
	int nsam;
{
//	double t;
	int i, *descl, *descr ;
	void parens( struct node *ptree, int *descl, int *descr, int noden );

	descl = (int *)malloc( (unsigned)(2*nsam-1)*sizeof( int) );
	descr = (int *)malloc( (unsigned)(2*nsam-1)*sizeof( int) );
	for( i=0; i<2*nsam-1; i++) descl[i] = descr[i] = -1 ;
	for( i = 0; i< 2*nsam-2; i++){
	  if( descl[ (ptree+i)->abv ] == -1 ) descl[(ptree+i)->abv] = i ;
	  else descr[ (ptree+i)->abv] = i ;
	 }
	parens( ptree, descl, descr, 2*nsam-2);
	free( descl ) ;
	free( descr ) ;
}

	void
parens( struct node *ptree, int *descl, int *descr,  int noden)
{
	double time ;

	if( descl[noden] == -1 ) {
	printf("%d:%5.3lf", noden+1, (ptree+ ((ptree+noden)->abv))->time );
	}
	else{
	printf("(");
	parens( ptree, descl,descr, descl[noden] ) ;
	printf(",");
	parens(ptree, descl, descr, descr[noden] ) ;
	if( (ptree+noden)->abv == 0 ) printf(");\n");
	else {
	  time = (ptree + (ptree+noden)->abv )->time - (ptree+noden)->time ;
	  printf(")%d:%5.3lf", noden+1, time );
	  }
		}
}

**  pickb : returns a random branch from the tree. The probability of picking
              a particular branch is proportional to its duration. tt is total
	      time in tree.   ***

	int
pickb(nsam, ptree, tt)
	int nsam;
	struct node *ptree;
	double tt;
{
	double x, y, ran1();
	int i;

	x = ran1()*tt;
	for( i=0, y=0; i < 2*nsam-2 ; i++) {
		y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
		if( y >= x ) return( i ) ;
		}
	return( 2*nsam - 3  );   changed 4 Feb 2010
}

	int
pickbmf(nsam, mfreq, ptree, tt )
	int nsam, mfreq;
	struct node *ptree;
	double tt;
{
	double x, y, ran1();
	int i, lastbranch = 0 ;

	x = ran1()*tt;
	for( i=0, y=0; i < 2*nsam-2 ; i++) {
	  if( ( (ptree+i)->ndes >= mfreq )  && ( (ptree+i)->ndes <= nsam-mfreq) ){
		y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
		lastbranch = i ;     changed 4 Feb 2010
	  }
	  if( y >= x ) return( i ) ;
	}
	return( lastbranch );     changed 4 Feb 2010
}

***  tdesn : returns 1 if tip is a descendant of node in *ptree, otherwise 0. *

	int
tdesn(struct node *ptree, int tip, int node )
{
	int k;

	for( k= tip ; k < node ; k = (ptree+k)->abv ) ;
	if( k==node ) return(1);
	else return(0);
}
*/


/* pick2()  */

	int
pick2(int n, int *i, int *j)
{
	double ran1();

	*i = n * ran1() ;
	while( ( *j = n * ran1() ) == *i )
		;
	return(0) ;
}

/**** ordran.c  ***/

/*
	void
ordran(int n,double pbuf[])
{
	ranvec(n,pbuf);
	order(n,pbuf);
	return;
}


	void
mnmial(int n, int nclass, double p[], int rv[])
{
	double ran1();
	double x, s;
	int i, j;

	for(i=0; i<nclass; i++) rv[i]=0;
	for(i=0; i<n ; i++) {
	   x = ran1();
	   j=0;
	   s = p[0];
	   while( (x > s) && ( j<(nclass-1) ) )  s += p[++j];
	   rv[j]++;
	   }
	return;
}

        void
order(int n,double pbuf[])
{
        int gap, i, j;
        double temp;

        for( gap= n/2; gap>0; gap /= 2)
           for( i=gap; i<n; i++)
                for( j=i-gap; j>=0 && pbuf[j]>pbuf[j+gap]; j -=gap) {
                   temp = pbuf[j];
                   pbuf[j] = pbuf[j+gap];
                   pbuf[j+gap] = temp;
                   }
    return;
}


	void
ranvec(int n,double pbuf[])
{
	int i;
	double ran1();

	for(i=0; i<n; i++)
		pbuf[i] = ran1();

	return;
}



	int
poisso(double u)
{
	double  cump, ru, ran1(), p, gasdev(double, double) ;
	int i=1;

	if( u > 30. ){
	    i =  (int)(0.5 + gasdev(u,u)) ;
	    if( i < 0 ) return( 0 ) ;
	    else return( i ) ;
	  }
	 
	ru = ran1();
	p = exp(-u);
	if( ru < p) return(0);
	cump = p;
	
	while( ru > ( cump += (p *= u/i ) ) )
		i++;
	return(i);
}


 a slight modification of crecipes version

double gasdev(m,v)
	double m, v;
{
	static int iset=0;
	static float gset;
	float fac,r,v1,v2;
	double ran1();

	if  (iset == 0) {
		do {
			v1=2.0*ran1()-1.0;
			v2=2.0*ran1()-1.0;
			r=v1*v1+v2*v2;
		} while (r >= 1.0);
		fac=sqrt(-2.0*log(r)/r);
		gset= v1*fac;
		iset=1;
		return( m + sqrt(v)*v2*fac);
	} else {
		iset=0;
		return( m + sqrt(v)*gset ) ;
	}
}
*/


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

#ifndef MS_NEW_H_
#define MS_NEW_H_


struct devent {
	double time;
	int popi;
	int popj;
	double paramv;
	double **mat ;
	char detype ;
	struct devent *nextde;
	} ;

struct c_params {
	int npop;
	int nsam;
	int *config;
	double **mig_mat;
	double r;
	int nsites;
	double f;
	double track_len;
	double *size;
	double *alphag;
	struct devent *deventlist ;
	} ;

struct m_params {
	 double theta;
	int segsitesin;
	int treeflag;
	int timeflag;
	int mfreq;
	 } ;

struct params { 
	struct c_params cp;
	struct m_params mp;
	int commandlineseedflag ;
	int output_precision;
	};

struct node{
	int abv;
	int ndes;
	float time;
	float brLength;
};

struct seg{
	int beg;
	int end;
	int desc;
	};

struct chromo{
	int nseg;
	int pop;
	struct seg  *pseg;
	} ;

struct segl {
        int beg;
        int end;
        struct node *ptree;
        int next;
        };

/*KRT -- prototypes added*/
void ordran(int n, double pbuf[]);
void ranvec(int n, double pbuf[]);
void order(int n, double pbuf[]);

void biggerlist(int nsam,  char **list );
int poisso(double u);
void locate(int n,double beg, double len,double *ptr);
void mnmial(int n, int nclass, double p[], int rv[]);
void usage();
int tdesn(struct node *ptree, int tip, int node );
int pick2(int n, int *i, int *j);
int xover(int nsam,int ic, int is);
int links(int c);

extern int brClass, mutClass, foldBrClass, allBrClasses, sampledPopsSize, poisTableSize;
extern int *poisTableScaleOrder;
extern int **blockLengthsMat;

#endif /* MS_NEW_H_ */
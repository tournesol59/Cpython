#include <stdio.h>
#include <math.h>
#include <malloc.h>

typedef struct IDENT05_TFdiscC {
    int nc;
    int nd;
    int size;
    double Ts;
    double *coeffs;
    double *yinit;
    //double *yin;
    //double *yout;
} ident05_tfdiscc;

ident05_tfdiscc* create_struct_tfdiscc(int dnc, int dnd, int dsize, double fTs, double* fcoeffs, double* yinit) 
{
    ident05_tfdiscc* tfdisc_s1 = (ident05_tfdiscc*) malloc(sizeof(ident05_tfdiscc));
    tfdisc_s1->nc = dnc;
    tfdisc_s1->nd = dnd;
    tfdisc_s1->size = dsize;
    tfdisc_s1->Ts = fTs;
    tfdisc_s1->coeffs = (double*) malloc(4*sizeof(double));
    tfdisc_s1->yinit =  (double*) malloc(2*sizeof(double));
}

double advanceintime(ident05_tfdiscc* tfdisc_sarg, int k) 
{
   int nc = tfdisc_sarg->nc;
   int nd = tfdisc_sarg->nd;
   double Ts = tfdisc_sarg->Ts;
   double ysum = 0;
   for (int i=0; i < nc; i++) 
   {
       ysum += sin((k-i)*Ts)*tfdisc_sarg->coeffs[i];
   }
   for (int i=0; i < nd; i++) 
   {
   // the first elmt in yinit is the most recent
       ysum += tfdisc_sarg->yinit[i]*tfdisc_sarg->coeffs[nc+i];    
   }
   return ysum;
}


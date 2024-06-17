#include <stdio.h>
#include <math.h>
#include "lapack.h"

int call_dgesv(int n, double *A, double *B) {
   lapack_int N=n;
   lapack_int RHS=n;
   lapack_int LDA=n;
   lapack_int * IPIV;
   lapack_int LDB=n;
   lapack_int INFO;
   int i;
 
   IPIV=(lapack_int*) malloc(n*sizeof(lapack_int));

//   INFO=LAPACKE_dgesv(LAPACK_COL_MAJOR,N,RHS,A,LDA, IPIV, B,LDB);
   INFO=lapack_dgesv(N,RHS,A,LDA,IPIV,B,LDB);
   if (INFO==1) {
     printf("lapacke error: dgesv has not found a solution");
     return -1;
   }
   else {
     return 0;
   }
}

int readvalues(double *a, double *b)  {

   FILE *fdesc;
//   int n=3;
//   int m=2;
//   int p=2;
   char **string=malloc(sizeof(char*) * 9); // allocate an array of 9 pointers  
   if (string==0) {
      fprintf(stderr,"Memory allocation failed");
      return -1;
   }
   int i = 0;
   //matrix a

   fdesc = fopen("matrices.txt", "r");
   if (fdesc == 0) {
      fprintf(stderr, "Could not open matrices.txt for reading");
      return -1;
   }
   for (i=0; i<9; i++) {
      string[i] = malloc(sizeof(char)*20);
      if (string[i]==0) {
         fprintf(stderr, "Memory allocation failed");
	 return -1;
      }
      fgets(string[i], 20, fdesc);
      a[i] = atof(string[i]); // let s try it
   }

//   fscanf("%lf %lf %lf", a[0], &a[1], &a[2]);
//   scanf("%lf %lf %lf", &a[3], &a[4], &a[5]);
//   fscanf("%lf %lf %lf", &a[6], &a[7], &a[8]);

   //matrix b
   for (i=9; i<12; i++) {
      string[i] = malloc(sizeof(char)*20);
      if (string[i]==0) {
         fprintf(stderr, "Memory allocation failed");
	 return -1;
      }
      fgets(string[i], 20, fdesc);
      b[i-9] = atof(string[i]); // let s try it
   }
//   fscanf("%lf %lf %lf", &b[0], &b[1], &b[2]);

   fclose(fdesc);   
   return 0;
}

int main(int argc, char *argv[]) {

   const int n=3;
   double a[9]={0.,0.,0.,0.,0.,0.,0.,0.,0.};   
   double b[3]={0.,0.,0.};
   double calloc[3]={0.,0.,0.};
   int i, sts;

   sts = readvalues(a, b);
   if (sts<0) {
        return -1;
   }
   sts = call_dgesv(n,a,b);
   if (sts <0) {
	return -1;
   }
   for (i=0; i<3; i++) {
      calloc[i]=b[i];
   }
//   printf("%lh %lh %lh\n", calloc[0], calloc[1], calloc[2]);

   printf("%lh \n", calloc[0]);
   
   return 0;
}

#define PY_SSIZE_T_CLEAN
#include "Python.h"
//essayer windows aussi
#include <stdio.h>
#include <math.h>
#include "lapacke.h"  // prefer set INCLUDE DIR in setup.py for portability

static PyObject *MatmultError = NULL;
static PyObject * matmult_system(PyObject *self, PyObject *args);


int call_dgesv(int n, double *A, double *B) {
   lapack_int N=n;
   lapack_int RHS=n;
   lapack_int LDA=n;
   lapack_int * IPIV;
   lapack_int LDB=n;
   lapack_int INFO;
   int i;
   for (i=0; i<n; i++) {
      IPIV[i]=0;
   }
   INFO=LAPACKE_dgesv(LAPACK_COL_MAJOR,N,RHS,A,LDA,IPIV,B,LDB);
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

/*
 * int savevalues(double *c) {
   FILE *fdesc;
//  int n=3;
//  fdesc = fopen("product.txt", "w");
//  fprintf("%lf %lf %lf", c[0], c[1], c[2]);
//  fprintf("%lf %lf", c[2], c[3]);
//  fprintf("%lf %lf", c[4], c[5]);
//  close(fdesc);
  return 0;
}
*
*/

int matmult(int n, int m, int p, double *a, double *b, double *c) {
   // matrices a,b,c are list of doubles
   // row_major_convention is assumed
   // dimensions are fixed: n=3, m=2, p=2 in waiting to know how to handle pointers
   int i,j,k;
   double s;

   for (k=0; k<p; k++) {
      for (i=0; i<n; i++) {
          s=0.;
	  for (j=0; j<m; j++) {
              s+= a[i*m+j]*b[j*p+k];
	  }
	  c[i*p+k]=s;
      }
   }
   if (c!=NULL) {
	return -1;
   }
   else {	
        return 0;
   }
}

static PyMethodDef MatmultMethods[] = {
   {"system", matmult_system, METH_VARARGS,
	    "numeric matrix multiplication"},
   {NULL,NULL,0,NULL}
};

static struct PyModuleDef matmultmodule = {
   PyModuleDef_HEAD_INIT,
   "matmult", // name of the module
    NULL, // module documentation, may be NULL (matmult_doc)
   -1,        // size of per-interpreter state of the module
              // or -1 if the module keeps state in global variables
   MatmultMethods  // defined above
};


PyMODINIT_FUNC PyInit_matmult(void) {
   PyObject *m;

   m = PyModule_Create(&matmultmodule);
   if (m==NULL) 
      return NULL;
   MatmultError = PyErr_NewException("matmult.error", NULL, NULL);
   Py_XINCREF(MatmultError);
   if (PyModule_AddObject(m, "error", MatmultError) < 0) {
       Py_XDECREF(MatmultError);
       Py_CLEAR(MatmultError);
       Py_DECREF(m);
       return NULL;
   }
   return m;
}


static PyObject * matmult_system(PyObject *self, PyObject *args)
{
   const int n=3;
   const int m=2;
   const int p=2;
//   PyObject *a;
   double a[9]={0.,0.,0.,0.,0.,0.,0.,0.,0.};   
//   PyObject *b;
   double b[3]={0.,0.,0.};
   double calloc[3]={0.,0.,0.};
   PyObject *arg2;
   int i, sts;

//   if (!PyArg_ParseTuple(args, "(dddddd)(dddd)", &a, &b))
//      return NULL;
//   if (!PyArg_ParseTuple(args, "(ddddddddd)(ddd)", &a, &b))
//      return NULL;
   sts = readvalues(a, b);
   if (sts<0) {
        PyErr_SetString(MatmultError, "numeric mult inputread failed");
        return NULL;
   }
   sts = call_dgesv(n,a,b);
   if (sts <0) {
	   PyErr_SetString(MatmultError, "numeric mult failed");
	   return NULL;
   }
   for (i=0; i<3; i++) {
       calloc[i]=b[i];
   }
   printf("%lf %lf %lf", calloc[0], calloc[1], calloc[2]);
   arg2 = Py_BuildValue("(ddd)", calloc); 
   Py_DECREF(a);
   Py_DECREF(b);
   return PyLong_FromLong(sts);
}



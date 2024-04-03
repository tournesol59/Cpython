#include "Python.h";

#include <math.h>

int matmult(int n, int m, int p, double *a, double *b, double *c) {
   // matrices a,b,c are list of doubles
   // row_major_convention is assumed
   // dimensions are fixed: n=3, m=2, k=2 in waiting to know how to handle pointers
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


static PyObject * matmult_system(PyObject *self, PyObject *args)
{
   const int n=3;
   const int m=2;
   const int p=2;
   PyObject *a;
   PyObject *b;
   PyObject *c;
   double calloc[6];
   
   int sts;

   if (!PyArg_ParseTuple(args, "(dddddd)(dddd)", &a, &b))
      return NULL;
   arg2 = Py_BuildValue("(dddddd)", calloc); 
   sts = matmult(n,m,p,a,b,arg2);
   if (sts <0) {
	   PyErr_SetString(MatmultError, "numeric mult failed");
	   return NULL;
   }
   Py_DECREF(a);
   Py_DECREF(b);
   return arg2;
}

static PyMethodDef MatmultMethods[] = {
   {"system", matmult_system, METH_VARARGS,
	    "numeric matrix multiplication"},
   {NULL,NULL,0,NULL}
};

static struct PyModuleDef matmultmodule = {
   PyModuleDef_HEAD_INIT,
   "matmult", // name of the module
   matmult_doc, // module documentation, may be NULL
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


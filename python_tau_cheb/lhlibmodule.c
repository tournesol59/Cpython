#include "python.h"
#include "math.h"
#include <stdio.h>
// 
static PyMethodDef LhLibMethods[] = {
    {"system",  lu_parser_wrapper, METH_VARARGS,
     "Execute a shell command."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


static struct PyModuleDef lhlibmodule = {
    PyModuleDef_HEAD_INIT,
    "Lhlib",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    LhlibMethods
};


PyMODINIT_FUNC PyInit_lhlib(void)
{
    return PyModule_Create(&lhlibmodule);
}

/* below the scientific part of our test module */
// commencer par une version sans prise en compte des pivots nuls
int inverse_lu_pivot(double *A, double *L, double *U, int m) {
   int ii, ij;
   // only for test: multtwo:
   for (ii=0; ii<arg1; ii++) {
      for (ij=0; ij<arg1; ij++) {
          U[ii*arg1+ij]=A[ii*arg1+ij];
      }
   }

   return 0; 
}


static PyObject *lu_parser_wrapper(PyObject *module, PyObject *args) 
{
   PyObject *ret = NULL;
   PyObject *pyStr = NULL;
   int arg1, arg2;
   FILE *inftab, *outfres; 
   double A[100];
   double L[100];
   double U[100];
   int ii,ij,res;
   PyObject *pyStr = NULL;
   arg2 = 8; /* optional arg Default value, when a max no lines for memalloc will be impl. */
   if (! PyArg_ParseTuple(args, "Si|i", &pyStr, &arg1, &arg2) {
      goto except;
   }
   /* your code here */
   inftab = fopen((const char*) pyStr, "r");
   for (ii=0; ii<arg1; ii++) {
      for (ij=0; ij<arg1; ij++) {
         fscanf(inftab, "%f", A[ii*arg1+ij]);
	 printf("%f ", A[ii*arg1+ij]);
      }
      fscanf(inftab, "\n");
      printf("%n");
   }
   fclose(inftab);

   res = inverse_lu_pivot(A, L, U, arg1);
   outfres = fopen("result", "w");
   for (ii=0; ii<arg1; ii++) {
      for (ij=0; ij<arg1; ij++) {
	 printf("%f ", L[ii*arg1+ij]);
         fprintf(outfres, "%f", L[ii*arg1+ij]);
      }
      fprintf(outfres, "\n");
   }
   
   for (ii=0; ii<arg1; ii++) {
      for (ij=0; ij<arg1; ij++) {
	 printf("%f ", U[ii*arg1+ij]);
         fprintf(outfres, "%f", U[ii*arg1+ij]);
      }
      fprintf(outfres, "\n");
   }

   fclose(outfres);

   /* end of code (provisory) */
   
   Py_INCREF(Py_None);
   ret=Py_None;
   assert(! PyErr_Occured());
   assert(ret);
   goto finally;
except:
   Py_XDECREF(ret);
   ret=NULL;
finally:
   return ret;

}

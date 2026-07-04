#include "Python.h"

static PyObject *basicsumstd_system(PyObject *self, PyObject *args)
{
   PyObject *item, *list, *squarelist;
   int i, n;
   PyObject *index;
   double xval, yval, meanval, stdval;
   PyObject *stdret;

   n=5; // defaut
   meanval=0.0;
   stdval=0.0;
   for (i=0; i <n; i++) {
      if (!PyAr_ParseTuple(args, "ddddd", xval))
         return -1;
      index = PyLong_FromSize_t(i);
      if (!index)
         return -1;
      item = PyFloat_FromDouble(xval);
      if (PyObject_SetItem(list, index, item) {
         PyDECREF(index);
         return -1;
      }
      meanval = meanval + xval;
      PyDECREF(index);
   }
   meanval = meanval / n;
// now calculating the list of squares
   for (i=0; i <n; i++) {
      PyObject *index = PyLong_FromSize_t(i);
      if (!index)
         return -1;
      if (PyObject_GetItem(list, index, item)) {
         PyDECREF(item);
         return -1;
      }
      xval = PyFloat_AsDouble(item, "d", xval);
      yval = (xval - meanval)*(xval - meanval);
      item = PyFloat_FromDouble(yval);
      if (PyObject_SetItem(squarelist, index, item) {
         PyDECREF(index);
	 PyDECREF(item);
         return -1;
      }
      stdval = stdval + yval 
      PyDECREF(index);
      PyDECREF(item);
   }
   stdval = stdval / n;
   stdret = PyFloat_FromFloat(stdval);
   
   return stdret;
}

static PyMethodDef BasicsumstdMethods[] = {
    {"system",  basicsumstd_system, METH_VARARGS,
     "Calculated standard deviation."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef basicsumstdmodule = {
    PyModuleDef_HEAD_INIT,
    "Basicsumstd",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    BasicsumstdMethods
};

PyMODINIT_FUNC PyInit_basicsumstd(void)
{
    return PyModule_Create(&basicsumstdmodule);
}



#include "Python.h"

PyObject *spam(void) {
   PyObject_Print(PyUnicode_FromFormat("Spam!"), stdout, 0);
   Py_RETURN_NONE;
}

PyMethodDef methods[]= {
  {"spam",
   (PyCFunction)spam,
   METH_NOARGS,
   "Say spam."},
   {NULL,NULL,0,NULL}};

PyModuleDef module = {
  PyModuleDef_HEAD_INIT,
  "spam",
  "Spam Module.",
  -1,
  methods};

PyMODINIT_FUNC PyInit_spam(void)
{
  return PyModule_Create(spammodule);
}


#define PY_SSIZE_T_CLEAN
#include "Python.h"
//essayer windows aussi

//first part: implementation of classes:
MATCLASS :: MATCLASS(int n, const char *matFile, const char *saveFile) :
   dim(n)
{
   strcpy(mat_file, matFile);
   strcpy(save_file, saveFile);
}

MATCLASS :: MATCLASS (const MATCLASS &M) 
{
   dim = M.dim;
   strcpy(mat_file, (const char*) M.mat_file);
   strcpy(save_file, (const char*) M.save_file);
}

MATCLASS :: ~MATCLASS() {}

const double * MATCLASS :: getC() {
   double tmp[3] = new double[3];
   for (int i=0; i<3; i++) {
      tmp[i]=C[i];
   }
   return tmp;
}

void MATCLASS :: setABC(const double *a, const double *b, const double *c) {
   for (int i=0; i<3; i++) {
      A[i]=a[i];
      B[i]=b[i];
      C[i]=c[i];
   }
}

int MATCLASS :: wrapper_dgesv() {
    for (int i=0; i<3; i++) {
       C[i]=B[i]; // because C will be erased and replaced by A*C
    }
    int sts = call_dgesv(dim, A, C);
    return sts;
}	
//second part: Cpp Extension module methods:
   // defining own types
/*
 * 
typedef struct {
   PyObject_HEAD
   PyObject *matfile;
   PyObject *savefile;
   int dim;
} CustomMatObject;
   // for defining own type, at a minimum we need a deallocation method 
static void CustomMat_dealloc(CustomMatObject *self)
{
   Py_XDECREF(self->matfile);
   Py_XDECREF(self->savefile);
   Py_TYPE(self)->tp_free((PyObject *) self);
}

  // custom new() method
static PyObject * Custom_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
   CustomMatObject *self;
   self = (CustomMatObject *) type->tp_alloc(type, 0);
   if (self != NULL) {
      self->matfile = PyUnicode_FromString(""):
      if (self->matfile == NULL) {
         Py_DECREF(self);
	 return NULL;
      }
      self->savefile = PyUnicode_FromString(""):
      if (self->savefile == NULL) {
         Py_DECREF(self);
	 return NULL;
      }
      self->dim=0;
   }
   return (PyObject *) self;
}
// we also define an initialization function
static int CustomMat_init(CustomMatObject *self, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"first", "last", "number", NULL};
    PyObject *matfile = NULL, *savefile = NULL, *tmp;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|00i", kwlist,
			   &matfile, &savefile,
			   &self->dim))
	    return -1;
    if (matfile) {
       tmp = self->matfile;
       Py_INCREF(matfile);
       self->matfile = matfile;
       Py_XDECREF(tmp);
    }
    if (savefile) {
       tmp = self->savefile;
       Py_INCREF(savefile);
       self->savefile = savefile;
       Py_XDECREF(tmp);
    }
    return 0;
}
// We might be tempted, for example, to assign thematfile member by INCRED self->matfile then DECREF matfile arg But this would be risky .. our type does not restrict the type of the matfile member, so itcould be any kind of object
static PyMemberDef CustomMat_members[] = {
   {"matfile", Py_T_OBJECT_EX, offsetof(CustomMatObject, matfile), 0, "matfile"},
    {"savefile", Py_T_OBJECT_EX, offsetof(CustomMatObject, savefile), 0, "savefile"},
    {"dim", Py_T_OBJECT_EX, offsetof(CustomMatObject, dim), 0, "dim"},
    {NULL}
};
// we define a single method, Custom.name() that outputs the objects name as the concatenation of the matfile and savefile names
static PyObject * CustomMat_name(CustomMatObject *self, PyObject *Py_UNUSED(ignored)) {
   if (self->matfile==NULL) {
      PyErr_SetString(PyExc_AttributeError, "matfile");
      return NULL;
   }
   if (self->savefile==NULL) {
      PyErr_SetString(PyExc_AttributeError, "savefile");
      return NULL;
   }
   return PyUnicode_FromFormat("%S %S", self->matfile, self->savefile);
}
// now that we defined the table: TBD

// table a completer
static PyTypeObject CustomType = {
//   .ob_base = PyVarObject_HEAD_INIT(NULL, 0),
//   .tp_name = "custom2.CustomMat",
//   .tp_doc = PyDoc_STR("Custom objects"),
//   .tp_basicsize = 0,
//   .tp_itemsize = 0,
//   .tp_flags = Py_TPFLAGS-DEFAULT | Py_TPFLAGS_BASETYPE,
//   .tp_new = CustoMat_new,
//   .tp_init = (initproc) CustomMat_init,
//   .tp_methods = CustomMat_methods
//   .tp_dealloc = (destructor) CustomMat_dealloc,
//   .tp_members = CustomMat_members,
//   .tp_methods = Custom_methods,
}

static PyMethodDef CustomMat_methods[] = {
   {"name", (PyCFunction) CustomMat_name, METH_NOARGS,
	    "numeric matrix custom type"},
   {NULL} //sentinel
};

static PyModuleDef customatmodule = {
   .m_base = PyModuleDef_HEAD_INIT,
   .m_name = "custom2",
   .m_doc = "Example module that creates an extension type.",
   .m_size = -1,
};
  *
  */

static PyObject *MatmultError = NULL;
static PyObject * matclass_system(PyObject *self, PyObject *args);

// precise that the function matclass_system can be added with MTH_VARARGS and METH_KEYWORDS
static PyMethodDef MatmultMethods[] = {
   {"system", matclass_system, METH_VARARGS | METH_KEYWORDS,
	    "numeric matrix multiplication"},
   {NULL,NULL,0,NULL}
};

static struct PyModuleDef matclassmodule = {
   PyModuleDef_HEAD_INIT,
   "matclass", // name of the module
    NULL, // module documentation, may be NULL (matclass_doc)
   -1,        // size of per-interpreter state of the module
              // or -1 if the module keeps state in global variables
   MatmultMethods  // defined above
};


PyMODINIT_FUNC PyInit_matclass(void) {
   PyObject *m;
i
 // PyObject *m2;
 // if (PyType_Ready(&CustomType) < 0)
 //    return NULL;
 //
 //    m2 = PyModule_Create(&custommatmodule);
 //    if (m2==NULL)
 //       return NULL;
 //
 //    if (PyModule_AddObject(m, "CustomMat", (PyObject *) &CustomType) < 0) {
 //       Py_DECREF(m2);
 //       return NULL;
 //    }

   m = PyModule_Create(&matclassmodule);
   if (m==NULL) 
      return NULL;
   MatmultError = PyErr_NewException("matclass.error", NULL, NULL);
   Py_XINCREF(MatmultError);
   if (PyModule_AddObject(m, "error", MatmultError) < 0) {
       Py_XDECREF(MatmultError);
       Py_CLEAR(MatmultError);
       Py_DECREF(m);
       return NULL;
   }
   return m;
}


static PyObject * matclass_system(PyObject *self, PyObject *args)
{
   PyObject *pyStr1 = NULL;
   PyObject *pyStr2 = NULL;
   int arg3;
   static char *kwlist[] = {
	   "matfile",
	   "savefile",
	   "dim",
	   NULL
   };
   PyObject *ret;

   char mat_file[8];
   char save_file[8];
//   const int n=3;
//   const int m=2;
//   const int p=2;
//   PyObject *a;
   double a[9]={0.,0.,0.,0.,0.,0.,0.,0.,0.};   
//   PyObject *b;
   double b[3]={0.,0.,0.};
   double calloc[3]={0.,0.,0.};
   int i, sts;

    if (! PyArg_ParseTupleAndKeywords(args, kwargs, "SSi", 
			   kwlist, pyStr1, &pyStr2, &arg3)) {
       goto except;
    }
   sts = strcpy(mat_file, (const char*) pyStr1);
   if (sts<0) {
        PyErr_SetString(MatmultError, "First String mult inputread failed");
        return NULL;
   }   
   sts = strcpy(save_file, (const char*) yStr2);
   if (sts<0) {
        PyErr_SetString(MatmultError, "Second String mult inputread failed");
        return NULL;
   }   
   sts = readvalues(mat_file,a, b);
   if (sts<0) {
        PyErr_SetString(MatmultError, "numeric mult inputread failed");
        return NULL;
   }
   matClassInst = new MATCLASS(arg3, mat_file, save_file);
//   if (!PyArg_ParseTuple(args, "(dddddd)(dddd)", &a, &b))
//      return NULL;
//   if (!PyArg_ParseTuple(args, "(ddddddddd)(ddd)", &a, &b))
//      return NULL;   matClassInst.setABC(a,b,calloc);

   // THE CORE OF THIS MODULE:
   sts = matClassInst.wrapper_dgesv(&colloc);  // operate on doubles inside the class
   if (sts <0) {
	   PyErr_SetString(MatmultError, "numeric mult failed");
	   return NULL;
   }
   // copy the result, though we could access it direct since it is public
   calloc = matClassInst.getC();
//   printf("%lf %lf %lf", calloc[0], calloc[1], calloc[2]);
   sts = savevalues((const char*) save_file, calloc);
   if (sts<0) {
        PyErr_SetString(MatmultError, "numeric mult outprint failed");
        return NULL;
   }

   ret = Py_BuildValue("(ddd)", calloc); 
   assert(! PyErr_Occurred());
   assert(ret);
   goto finally;

   except:
      Py_XDECREF(ret);
      ret = NULL;
//   return PyLong_FromLong(sts);
   finally:
      return ret;
}



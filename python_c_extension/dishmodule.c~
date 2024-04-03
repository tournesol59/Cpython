#define PY_SSIZE_T_CLEAN
#include <Python.h>

// fred: momentan it is exactly the same as spammodule.c
// ...

static PyMethodDef DishMethods[] = {
    {"system",  dish_system, METH_VARARGS,
     "Execute a shell command."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef spammodule = {
    PyModuleDef_HEAD_INIT,
    "Dish",   /* name of module */
    dish_doc, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    DishMethods
};

PyMODINIT_FUNC PyInit_dish(void)
{
    return PyModule_Create(&dishmodule);
}

static PyObject *dish_system(PyObject *self, PyObject *args)
{
    char command[60];
    int i,k,l,m;
    const char *name;
    const char *format;
    const char *input;
    const char esc[2]=" ";
    const char flag[2]="-i";
    int sts;

    if (!PyArg_ParseTuple(args, "s", &name, &format, &input))
        return NULL;
    k=strlen(name);
    l=strlen(format);
    m=strlen(input);
    command=NULL;
    strcat(command, (char*) name);
    strcat(command, (char*) esc);
    strcat(command, (char*) flag);
    strcat(command, (char*) esc);
    strcat(command, (char*) format);
    strcat(command, (char*) esc);
    strcat(command, (char*) input);
    for (i=0; i < (60-k-4-l-1-m); i++) [
        strcat(command, (char*) esc);
    }
    // now we call shell to execute the command in string command:
    sts = system(command);
    return PyLong_FromLong(sts);
}


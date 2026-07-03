#define PY_SSIZE_T_CLEAN
#include <Python.h>

// fred: momentan it is exactly the same as spammodule.c
// ...
static PyObject *dish_system(PyObject *self, PyObject *args)
{
    char command[80];
    int i,k,l,m;
    const char name[23]="echo 'HelloWorld' | sed";
    const char esc[2]=" ";
    const char *format;
    const char *input;
    const char *flag;
    //char format[15];
    //char input[20];
    //const char flag[12]="--expression=";
    int sts;

    printf("Parsing the argumentsi\n");

    if (!PyArg_ParseTuple(args, "sss", &flag, &format, &input))
        return NULL;
    
    k=strlen(flag);
    l=strlen(format);
    m=strlen(input);
    printf("Beginning of cmd: %s\n", (char*) &name[0]);
    printf("Argument flag: %s \n", &flag[0]);
    printf("Argument format: %s \n", &format[0]);
    printf("Argument filename: %s \n", &input[0]);
    // concatenate in series into the "command" char
    // eg it will be similar to "sed  's/kong/hong' hong-kong

    strcpy(command, "");
    strcat(command, (char*) name);
    strcat(command, (char*) esc);
    strcat(command, (char*) flag);
//    strcat(command, (char*) esc);
    strcat(command, (char*) format);
    // try pipe for the moment, input file part commented
    //strcat(command, (char*) esc);
    //strcat(command, (char*) input);
    printf("completing the char string\n");
    // fills the rest characters
    //for (i=0; i < (60-3-1-k-1-l-1-m); i++) {
    // for (i=0; i < (80-21-1-k-l-5); i++) {
    for (i=0; i < (80-23-1-k-l-5); i++) {
        strcat(command, (char*) esc);
    }
    // now we call shell to execute the command in string command:
    printf("now we have the command: %s\n", command);
    sts = system(command);
//    sts = system("echo 'HelloWorld!' | sed --expression=s/Hello/hello/g"); 
//                 			 ^		^		^
//					23		k		l
    sts=0;
    return PyLong_FromLong(sts);
}

static PyMethodDef DishMethods[] = {
    {"system",  dish_system, METH_VARARGS,
     "Execute a shell command."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef dishmodule = {
    PyModuleDef_HEAD_INIT,
    "Dish",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    DishMethods
};

PyMODINIT_FUNC PyInit_dish(void)
{
    return PyModule_Create(&dishmodule);
}


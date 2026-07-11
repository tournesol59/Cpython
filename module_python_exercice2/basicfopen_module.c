#include <Python.h>
// programme en vue dun expose sur la transmission de parametres variables dans un tampon (const char*) depuis Python vers le C
// l utilisation comme pour la base de donnees: necessite un tuple [database, user, passwd]

#define PY_SSIZE_T_CLEAN
#include <Python.h>

// ...
// declaration

char **str_split(char * a_str, const char a_delim);

static PyObject *basicfopen_system(PyObject *self, PyObject *args)
{
/*
    char command[81];
    command[0]='\0'; 
    int i,k,l,m;
    //const char name[24]="echo 'HelloWorld' | sed";
    char name[24];
    name[0]='\0';
    const char name1[6]="echo '";
    const char name3[7]="' | sed";
    const char esc[2]=" ";
    const char *format;
    const char *input;
    const char *inputdeb="HelloWorld";
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
*/

  //input as const char*
    const char *input;
    const char *modus;
  // size of the input string
    size_t input_len;
    int nargs;

    if (!PyArg_ParseTuple(args, "ss", &input, &modus))
        return NULL;
    
    input_len=strlen(input);
    //

    const char a_delim=","; 
    //use of the basic library str_split function
    // return an array of pointer to char **
    char** input_split;
    char filename[50]={0};
    char username[50]={0};
  
    input_split = str_split((char*) input, a_delim);
    strcpy(filename, (const char *) input_split[0]);
    printf("filename: %s\n", filename);

    strcpy(username, (const char *) input_split[1]);
    printf("username: %s\n", username);

    PyObject* result = Py_BuildValue("sss", filename, username, "pass");

    return result; // try this, if error got, search for PyUnicode 
}

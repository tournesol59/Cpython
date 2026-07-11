// programme en vue dun expose sur la transmission de parametres variables dans un tampon (const char*) depuis Python vers le C
// l utilisation comme pour la base de donnees: necessite un tuple [database, user, passwd]

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <regex.h>
#define MAX_MATCHES 3
// ...
// declaration

char **str_split(char * a_str, const char a_delim);

static PyObject *basicfopenquery_system(PyObject *self, PyObject *args)
{
  /**** preliminary section for parsing the input(url)  ***/ 
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


    // use input string the regex POSIX library
    char filename2[50]={0};
    char username2[50]={0};
    const char *pattern="sqlite\.(.*)\.(.*)";  //"sqlite.sample.fred"; 
    regex_t *preg;
    regmatch_t matches[MAX_MATCHES];
    // compile the pattern 
    if (!regcomp(preg, pattern, REG_NEWLINE)) {
       fprintf(stderr, "regcomp error check the regex pattern");
    }
    //use the pattern
    if (regexec(preg, input, 3, matches, 0)) {
       // success
      memcpy(filename2, input+matches[1].rm_so, matches[1].rm_eo - matches[1].rm_so +1);
      memcpy(username2, input+matches[2].rm_so, matches[2].rm_eo - matches[2].rm_so +1);
    }
    else
      fprintf(stderr, "regexec error, either the input string is mismatched or the pattern wrongly explicited"); 


    //now opens the file, whose name was transmitted and read its first line:
    if (! openforread_fgets(filename)) {
       fprintf(stderr, "Error openforread_fgets: check the filename %s\n", filename);
    } 
    PyObject* result = Py_BuildValue("sss", filename, username, "pass");

    return result; // try this, if error got, search for PyUnicode 
}

int openforread(char *filename) {
   int is_ok = EXIT_FAILURE;
   FILE* fp = fopen(filename, "r");
   char buffer[1000];
   if (!fp)
   {
      fprintf(stderr, "File opening failed");
      return is_ok;
   }
   int i=0;
   int c; // Note: int, not char, required to handle EOF 
   while ((c = fgetc(fp)) != '\\0') {
      buffer[i]=c;
      i++;
   } 

   printf("first line: %s\n", buffer);
   
   if (ferror(fp))
      puts("I/O error when reading");
   else if (feof(fp))
   {
      puts("End of file reached successfully");
      is_ok = EXIT_SUCCESS;
   }
   fclose(fp);
   return is_ok;
}


int openforread_fgets(char *filename) {
   int is_ok = EXIT_FAILURE;
   FILE* fp = fopen(filename, "r");
   char buffer[15];
   if (!fp)
   {
      fprintf(stderr, "File opening failed");
      return is_ok;
   }
   int i=0;
   while (1) {
      if (!fgets(buffer, sizeof buffer, fp)) break; 
      printf("line %d : %s\n", i, buffer);
   }
   if (feof(fp)) {
      fprintf(stderr, "End of file reached");
      is_ok = EXIT_SUCCESS;
   } 
   fclose(fp);
   return is_ok;
} 

static PyMethodDef BasicfopenqueryMethods[] = {
    {"system",  basicfopenquery_system, METH_VARARGS,
     "Execute a shell command."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef basicfopenquerymodule = {
    PyModuleDef_HEAD_INIT,
    "basicfopenquery",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    BasicfopenqueryMethods
};

PyMODINIT_FUNC PyInit_basicfopenquery(void)
{
    return PyModule_Create(&basicfopenquerymodule);
}


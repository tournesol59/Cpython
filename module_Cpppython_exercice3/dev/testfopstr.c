#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "idresult.h"
#include "fopstrset.h"
// une petite structure de pile elementaire qui aide
// voir "idresult.c" et "libfopstrset.c"

int main(int argc, char **argv) {
   
   const char *filename="sampledb";
   char *result;
   //const char* url="sqlite.sample.alan";
   int i=0;
   char **result_t=malloc(5*sizeof(char*));
   char item[40];
   //first of all test of our implementation of queue result_set{ 
   result_setp resultLL_p;
   
   //then core of our prgm
   if (openforreadline_s((char*) filename, &resultLL_p)==EXIT_FAILURE) {
     printf("error opening of file sample");
   }

   while (resultLL_p != NULL) {
      printf("item %d: ",i);
      pop(resultLL_p);
      i++;
   }
 
   return 0;
}



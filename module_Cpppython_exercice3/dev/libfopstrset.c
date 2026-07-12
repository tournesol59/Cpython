/* sub programs and function to parse a text file separated with commas
 * and sanitize it
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "idresult.h"
#include <assert.h>

int validstrsep(char *buffer, const char* a_delim, int count1, char **result) {
   int is_ok=EXIT_FAILURE;
   char delim[2];
   int count;
   char *tmp;
   char *last_comma;

   //ssize_t strmax = sizeof(buffer);
   char* token;
   char* next_token;
   //char** result;
   int idx;

   //step 1, assign separator and count the number of separators
   strncpy(delim,a_delim,1);
   delim[1]=0; 
   tmp=buffer;
   count=0;
   while(*tmp) {
     if (!strncmp(tmp, (const char*) delim, 1)) {
        count++;
     }
     tmp++; 
   }
   // increase if a seperator is at the end
   count += last_comma < (buffer+strlen(buffer)-1);
   // add space for the terminating null string
   count++;
 
   if ((count1!=0) && (count1!=count)) {
     is_ok=EXIT_FAILURE;
     return is_ok; 
   } 
  
   //step 2 , in a while(token) get the token address and save its content:
   printf("Parsing the input string: '%s'\n", buffer);
   result=(char**) malloc((count+1)*sizeof(char*));
   idx=0;
   //token = strtok_s(buffer, &strmax, delim, &next_token);
   token = strtok(buffer, delim);
   while (token)
   {
     assert(idx<count-1); 
     *(result+idx)=strdup(token);
     //token=strtok_s(NULL, &strmax, delim, &next_token);
     token=strtok(NULL, delim);
     idx++;     
   }
   assert(idx==(count-1)); //count-1; 
   *(result+idx++)=0;

   is_ok=0;
   return count;
}

//int openforreadline(char *filename);
int openforreadline_s(char *filename, result_setp* resLL_p) {
   int is_ok=EXIT_FAILURE;
   char buffer[100];
   const char *a_delim=",";

   size_t count=0;
   char *tmp;
   char *last_comma;
   int irow=0;
   int idx=0;
  
   char** result;

   //step 1 opening the file
   FILE *fp = fopen(filename, "r");
   if (!fp) {
	fprintf(stderr, "file opening failed");
        return is_ok;
   }
   // step2: loop over sthe readine
   count=0; 
   while (fgets(buffer, sizeof buffer, fp) != NULL) {  
//      if (fgets(buffer, sizeof buffer, fp) == NULL) {  

      // step2-1 delegate the task of counting and extracting tokens inside separators
      //      to the subroutine validstrsep
      count=validstrsep(buffer, a_delim, count, result);

     // step 3: copy the pointed to char in a structure linked list:
      printf("name: %s, age: %d, body: %s\n", *(result),*(result+1),*(result+2)); //debug
      char name[40];
      strcpy(name, (const char*) *(result+0));
      int age = atol((const char*) *(result+1));
      char body[40];
      strcpy(body, (const char*) *(result+2));
      printf("name: %s, age: %d, body: %s\n", *(result),*(result+1),*(result+2)); //debug

      if (irow==0) {
         *resLL_p = create(name, age, body);
      }
      else {
         *resLL_p = add(name, age, body, *resLL_p);
      }
      printf("VERIF name: %s, age: %d, body: %s\n", (char*) (*resLL_p)->name, (*resLL_p)->age, (char*) (*resLL_p)->body); //debug
      irow++; 
   }
   fclose(fp); 
   is_ok=0;
   return is_ok;
}



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

typedef struct result_set{
    int id;
    char name[40];
}; 
//int openforreadline(char *filename);
int openforreadline(char *filename, char **result) {
   int is_ok=EXIT_FAILURE;
   char buffer[100];
   const char *a_delim=",";
   char delim[2];
   size_t count=0;
   char *tmp;
   char *last_comma;
   int i=0;
   int idx=0;
//   int j=0;
   strcpy(delim,a_delim);
   delim[1]=0;
  
   //step 1 opening the file
   FILE *fp = fopen(filename, "r");
   if (!fp) {
	fprintf(stderr, "file opening failed");
        return is_ok;
   }
   // step2: loop over sthe readine 
   while (fgets(buffer, sizeof buffer, fp) != NULL) {  

//      if (fgets(buffer, sizeof buffer, fp) == NULL) {  

      // step2-1 count how many element will be extracted
      count=0;
      tmp=buffer;
      while (*tmp) {
         if (!strncmp(tmp,a_delim,1)) { // this side: ==0 means there is a delimiter i
            count++;
            last_comma=tmp;
	 }
	 tmp++;
      }
      // increase if a seperator is at the end
      count += last_comma < (buffer+strlen(buffer)-1);
      // add space for the terminating null string
      count++;
      printf("number of token expected: %d\n", count);
     // step 2-2: first call to strtok, changed to (count+1) after tests 
      //result=malloc((count+1)*sizeof(char*)); 
      if (result) {
	 idx=0; 
         char *token = strtok(buffer, delim);
         while ((token) && (idx < (count-1))) {   
            assert(idx < count-1); //count-1;
	   *(result+idx++)=strdup(token); // see the doc
	 //step 2-3x:subsequent calls to strtok this timewith nul inputs
	    token=strtok(0, delim);
	 }
	 assert(idx==(count-1)); //count-1; 
         *(result+idx++)=0; 
      }
     i++; 
   }
   fclose(fp); 
   is_ok=0;
   return is_ok;
}

int openforreadline_separator(char *filename, char **result) {
   int is_ok=EXIT_FAILURE;
   char buffer[100];
   const char *a_delim=",";
   char delim[2];
   strcpy(delim,a_delim);
   delim[1]=0;
   int i=0;
   int idx=0;
   char *pstr; 
   //step 1 opening the file
   FILE *fp = fopen(filename, "r");
   if (!fp) {
	fprintf(stderr, "file opening failed");
        return is_ok;
   }
  while (i<1) {
   // step2: loop over sthe readine 
      if (fgets(buffer, sizeof buffer, fp) == NULL) {  
         is_ok = EXIT_FAILURE;
	 fprintf(stderr, "openforrreadline error fgets");
	 return is_ok; 
      }
      idx=0;
      pstr=strsep(&buffer, delim);
      strcpy(*(result+idx),(const char*) pstr);
      idx++; 
      while (result != NULL) {
         pstr=strsep(&buffer, delim);
         strcpy(*(result+idx),(const char*) pstr);
	 idx++;
      } 

      i++;
   }
   fclose(fp); 
   is_ok=0;
   return is_ok;
}


int main(int argc, char **argv) {
   
   const char *filename="sampledb";
   char *result;
   //const char* url="sqlite.sample.alan";
   int i;
   char **result_t=malloc(5*sizeof(char*));
   char *item[40];
   if (openforreadline((char*) filename,result_t)==EXIT_FAILURE) {
     printf("error opening of file sample");
   }
   printf("line: ");
   for (i=0; i<3; i++) {
     strcpy(item, (const char*) *(result_t+i));
     printf("%s, ", item); 
   }
   printf("\n");
   return 0;
}



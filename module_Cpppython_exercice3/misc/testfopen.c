#include <stdio.h>
#include <stdlib.h>
#include <regex.h>
#include <string.h>

#define MAX_MATCHES 5 
#define MAX_BUFFER 100

#define PRGM_CLASS_OPENFILES testfopen_

int handleregex(char *input, char *target);

int openforread(char *filename);

int openforreadbuf(char *filename);

void ReportException(char *subname, int err_code);

int main(int argc, char **argv) {
   
   const char *filename="sample";
   char *result;
   const char* url="sqlite.sample.alan";
   int i;

   if (openforreadbuf((char*) filename)==EXIT_FAILURE) {
     printf("error opening of file sample");
   }

   result = (char*) malloc(10*MAX_BUFFER*sizeof(char));
   if (handleregex((char*) url, result)) {
     for (i=0;i<2;i++) {
       printf("matched string i=%d: %s\n", i, result+i*100);
     }      
   }
   return 0;
}

int handleregex(char *input, char *target) {
  const char *pattern="sqlite\\.(.*)\\.(.*)"; 
  regex_t preg;
  //char *target;
  int i;
  int is_ok=1; 
  if (regcomp(&preg, pattern, REG_EXTENDED)) {
     printf("error compile regex");    
     is_ok=0; 
  } 
  int nmatch=preg.re_nsub;
  regmatch_t pmatch[MAX_MATCHES];
  if (!regexec(&preg, input, nmatch+1, pmatch, 0)) {
    for (i=1; i<nmatch+1; i++) {
       memcpy(&target[i*MAX_BUFFER], input+pmatch[i].rm_so, pmatch[i].rm_eo - pmatch[i].rm_so +1); 
       target[i*MAX_BUFFER+pmatch[i].rm_eo - pmatch[i].rm_so]='\0';
    } 
  }
  else {
    printf("error execution regex");
    is_ok=0;
  }
  return is_ok; 
}

int openforread(char *filename) {
   int is_ok = EXIT_FAILURE;
   FILE* fp = fopen(filename, "r");
   char buffer[1000];
   if (!fp)
   {
      printf("File opening failed");
      return is_ok;
   }
   int i=0;
   int c; // Note: int, not char, required to handle EOF 
   while ((c = fgetc(fp)) != '\n') {
      buffer[i]=c;
      i++;
   } 

   printf("first line: %s\n", buffer);
   
   if (ferror(fp))
      puts("I/O error when reading");
   else if (feof(fp))
   {
      printf("End of file reached successfully");
      is_ok = EXIT_SUCCESS;
   }
   fclose(fp);
   return is_ok;
}

/*** other implementation, based on fgets()  ***/

int openforreadbuf(char *filename) {
   int is_ok = EXIT_FAILURE;
   FILE* fp = fopen(filename, "r");
   char buffer[15];
   if (!fp)
   {
      printf("error File opening failed (fgets)");
      return is_ok;
   }
   int i=0;
   while (1) {
      if (!fgets(buffer, sizeof buffer, fp)) break; 
      printf("line %d : %s\n", i, buffer);
      i++;
   }
   if (feof(fp)) {
      printf("End of file reached");
      is_ok = EXIT_SUCCESS;
   } 
   fclose(fp);
   return is_ok;
} 

void ReportException(char *subname, int err_code) 
{
   char packgname[100]={0};
   strcat(packgname, "testfopen_"); //'PRGM_CLASS_OPENFILES'
   strcat(packgname, "openforreadbuf");
   if (!strncmp(subname, "openforreadbuf", 14)) 
   {
      printf("%s: %d the file may be corrupt or not in the directory", packgname, err_code);       
   }
}


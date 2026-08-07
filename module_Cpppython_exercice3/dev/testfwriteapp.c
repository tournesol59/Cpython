#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int openforreadandwrite(char* filename) {

// open file in append+ mode, as for the documentation, it should open in read access
// and enabled to write to end of the file, but output cannot be followed by input without
// an intervening call to fflush
   int is_ok = EXIT_FAILURE;
   char buffer[100];
   const char* a_delim=",";
   char delim[2];
   const char* a_delend="\n";
   char delend[2];
   const char* name="Sophie Bell";
   int age=30;
   char agestr[3];
   const char* speciality="Medicine";
 
   FILE *fs;
   
   fs = fopen(filename, "a+");
   // append operation:  
   // first step is to concatenate the words into a sentence, buffer
   strcpy(delim, a_delim);
   delim[1]=0;
   strcpy(delend, a_delend);
   delend[1]=0;
   
   strcat(buffer, name);
   strcat(buffer, (const char*) delim);
   sprintf(agestr, "%d", age);
   strcat(buffer, agestr);
   strcat(buffer, (const char*) delim);
   strcat(buffer, speciality);
   strcat(buffer, delim);

   strcat(buffer, delend);
   
   // second step is to append to a file
   /*
   if (fputs((const char*) buffer, fs) <0) 
   {
      printf("error cannot write appended string to the test file sampledb4w\n");
      is_ok = 1;
      return is_ok;  
   }
    // before read operation, we need to call fflush
    if (fflush(fs)!=0) 
    {
       printf("error cannot flush the internal buf of the test file sampledb4w\n");
       is_ok=2;
       return is_ok;
    }
    */
     // read operation
     
    rewind(fs); 
    char inbuf[100];
    while ((fgets(inbuf, sizeof(inbuf), fs)) != NULL) 
    {
       printf("%s \n", inbuf);   
    }
    if (feof(fs)) 
    {
       printf("end of file reached");   
    }
    
    // close the file
    fclose(fs);
    is_ok=0;
    return is_ok;
}

int main(int argc, char** argv) {
   const char* filesample = "sampledb4w";
   char filename[20];
   
   strcpy(filename, filesample);
   if (openforreadandwrite(filename)) 
   { 
      printf("error cannot open file for edit \n");   
   }
   
   return 0;
}
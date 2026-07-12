#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct result_set {
   int id;
   char *name;
};

int glb_iid = 0;

int openforsearch(FILE* fp, char *searchit, struct result_set *item) {
   int is_ok = EXIT_FAILURE;
//   FILE* fp = fopen(filename, "a+");
   char buffer[100];
   int iid;
   char iname[40];
   int i=0;
   while (1) {
      if (!fgets(buffer, sizeof buffer, fp)) break; 
      //printf("line %d : %s\n", i, buffer);
      // process that line;
      sscanf(buffer, "%d, %s", &iid, iname);
      if (strcmp(iname, (const char*) searchit)) {
         // found item that matches name
	 item->id = iid;
         strcpy(item->name, (const char*) iname);	
         //break; 
      }
      i++;
   }
   if (feof(fp)) {
      printf("End of file reached");
      is_ok = EXIT_SUCCESS;
   }
// set the max number of rows in the file as global
   glb_iid = i; 
  // fclose(fp);
   return is_ok;
} 

int openforappend(FILE* fp, struct result_set *item) {
   int is_ok = EXIT_FAILURE;
//   FILE* fp = fopen(filename, "a+");
   int iid=2;
   char iname[40];
   char buffer[100];
   strcpy(iname, "Fran All");
   sprintf(buffer, "%d, %s", iid, iname);
   if (fputs((const char*) buffer, fp)) {
      printf("error write to file");
      is_ok=2;
      return is_ok;
   }
   is_ok=0;
   return is_ok;

}

//int openforsearch(FILE* fp, char *searchit, struct result_set *item) {
//int openforappend(FILE* fp, struct result_set *item);
 
int main(int argc, char** argv) {

   const char *filename = "sample";
   const char *nameit = "Alan Tur";
   struct result_set result_it;

  FILE* fp = fopen(filename, "a+");
  if (!fp)
   {
      printf("error File opening failed (fopen)");
      return 1;
   }
 
  if (openforsearch(fp, (char*) nameit, &result_it)) {
      printf("error openforappend");
      exit(1); 
   }
   printf("result set: %d, %s\n", result_it.id, result_it.name);

   fclose(fp);
   return 0;
}



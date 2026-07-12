#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

// une petite structure de pile elementaire qui aide
// a resoudre la probleme de tableau de pointeur de pointeurs
typedef struct result_set RESULT_SET;
struct result_set {
    int id;
    char name[40];
    int age;
    char body[40];
    RESULT_SET *next;
};
typedef RESULT_SET* result_setp;
//create

RESULT_SET* create(char *name, int age, char *body) {
   RESULT_SET *ptr=(result_setp) malloc(sizeof(RESULT_SET));
   ptr->id=1;
   strcpy(ptr->name, (const char*) name);
   ptr->age=age;
   strcpy(ptr->body, (const char*) body);
   ptr->next=NULL;
   return ptr;
}
//add on top
result_setp add(char *name, int age, char *body, RESULT_SET *head) {
   result_setp ptr=create(name,age,body);
   //RESULT_SET rs=create(name,age,body);
   //result_setp ptr=&rs;
   ptr->next = head;
   return ptr;
}
//pop from top: display and remove
result_setp pop(RESULT_SET *ptr) {
   if (ptr==NULL) {
     exit(1);
   }
   else {
     result_setp pnext_c = ptr->next;
     printf("name: %s, age: %d, body: %s\n", ptr->name,ptr->age,ptr->body);
     free(ptr);
     return pnext_c;
   }
}

//int openforreadline(char *filename);
int openforreadline(char *filename, result_setp* resLL_p) {
   int is_ok=EXIT_FAILURE;
   char buffer[100];
   const char *a_delim=",";
   char delim[2];
   size_t count=0;
   char *tmp;
   char *last_comma;
   int irow=0;
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
         if (!strncmp(tmp,a_delim,1)) { // this side: ==0 means there is a delimiter
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
      char** result=malloc((count+1)*sizeof(char*)); 
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
      pstr=strsep((char**) &buffer, delim); // added a cast
      strcpy(*(result+idx),(const char*) pstr);
      idx++; 
      while (result != NULL) {
         pstr=strsep((char**) &buffer, delim); // added a cast
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
   int i=0;
   char **result_t=malloc(5*sizeof(char*));
   char item[40];
   //first of all test of our implementation of queue result_set{ 
   result_setp resultLL_p;
   
   //then core of our prgm
   if (openforreadline((char*) filename, &resultLL_p)==EXIT_FAILURE) {
     printf("error opening of file sample");
   }

   while (resultLL_p != NULL) {
      printf("item %d: ",i);
      pop(resultLL_p);
      i++;
   }
   /*
   printf("line: ");
   for (i=0; i<3; i++) {
     strcpy(item, (const char*) *(result_t+i));
     printf("%s, ", item); 
   }
   printf("\n");
   */
   return 0;
}



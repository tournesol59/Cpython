#include "idqueue_adv.h"

#define IDRESULT_EXIT_FAILURE 41

//length
int queue_length(result_setp head) {
   int count=0;
   result_setp tmp = head;
   while (tmp != NULL) {
      count++;
      tmp = tmp->next;
   }
   return count;
}

// init
result_setp queue_create(const char *name, int age, const char *body) {
   RESULT_SET* ptr = (RESULT_SET*) malloc(sizeof(RESULT_SET));
   if (ptr != NULL) {
      ptr->name = (char*) malloc(sizeof(name));
      ptr->body = (char*) malloc(sizeof(body));
      if ((ptr->name != NULL) && (ptr->body != NULL)) {
         strcpy(ptr->name, (const char*) name);
         strcpy(ptr->name, (const char*) body);   
         ptr->age = age;
         ptr->next=NULL;
      }
      else {
         return NULL;
      }
   }
   return ptr;
}

//add at end position
int queue_addlast(const char *name, int age, const char *body, RESULT_SET* head)
{
   int exit_ok=-1;   
   RESULT_SET *tmp=head;
   RESULT_SET *tmp_ins = queue_create(name, age, body);
   if (head == NULL) {
      head = tmp_ins;
   }
   else if (head->next == NULL) {
      head->next = tmp_ins;
   }
   else {
      while (tmp->next != NULL) {
         tmp = tmp->next;
      }
      tmp->next = tmp_ins;
   }
   
   exit_ok=0;
   return exit_ok;
}

//update match on name
int queue_updatematch(const char *matchname, int newage, char* newbody, RESULT_SET *head)
{
   int exit_ok=-1;
   RESULT_SET* tmp = head;
   if (tmp == NULL) {
      exit_ok=1;
      return exit_ok;
   }
   else {
      while ((tmp->next != NULL) || (strcmp(tmp->name, matchname) != 0)) {
         tmp=tmp->next;
      }
      if (tmp->next != NULL) {
         tmp->age = newage;
         strcpy(tmp->body, (const char*) newbody);
         exit_ok=0;
      }
      else if (!strcmp(tmp->name, matchname)) {
         // here we must do the same thing but we could not be sure w/o the test second time
         tmp->age = newage;
         strcpy(tmp->body, (const char*) newbody);
         exit_ok=0;
      }
      else {
         exit_ok=2;
      }
   }
   return exit_ok;
}

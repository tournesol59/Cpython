#include "idresult.h"
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

//search for a particular cell for cell->name ="matchname"
// and replace its body and age
result_setp replace(char* matchname, int newage, char* newbody, RESULT_SET *head) {
   RESULT_SET* tmp;
   RESULT_SET* tmp_prev;
   tmp=head;
   tmp_prev=head;
   int count;
   if (tmp==NULL) {
      printf("empty linked list");
      exit(1);
   } else {
      while ((tmp != NULL) && (strcmp(tmp->name, (const char*) matchname))) {
         tmp_prev = tmp;
         tmp = tmp->next;
         count++;
      }
      if (tmp !=NULL) {//matched item
         RESULT_SET* tmp_ins = create(tmp->name, newage, newbody);
         tmp_ins->next = tmp->next;
         //particular case: matched item on head, that why we used a counter
         if (count==0) {
            head=tmp_ins;
         } else {
            tmp_prev->next=tmp_ins;
         } 
      }
   }
   return head;
}
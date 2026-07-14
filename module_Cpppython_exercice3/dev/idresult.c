#include "idresult.h"

//list_length()
int list_length(result_setp head) {
   RESULT_SET *ptr=head;	
   int count=0;
   while (ptr != NULL) {
     ptr=ptr->next;
     count++;
   }
   return count;   
}
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
result_setp replace(char* matchname, int newage, char* newbody, int selected_id, RESULT_SET *head) {
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
      if (tmp != NULL) {//matched item
       if (tmp->id == selected_id) {
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
   }
   return head;
}

//search for a particular cell for cell->name ="matchname"
//and return another linked list with all matching cell for field cell->name
result_setp searchname(char* matchname, RESULT_SET *head) {
   RESULT_SET* tmp;
   RESULT_SET* tmp_prev;
   RESULT_SET* tmp_match;
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
      }
      if (tmp != NULL) {//matched item
        if (count==0) { 
	 tmp_match = create(tmp->name, tmp->age, tmp->body);
         tmp_match->id =tmp->id;
         count++;
        }
	else {
         tmp_match = add(tmp->name, tmp->age, tmp->body, tmp_match);
	 tmp_match->id =tmp->id;
         count++;
	}
      }
   }
   return tmp_match;      
}

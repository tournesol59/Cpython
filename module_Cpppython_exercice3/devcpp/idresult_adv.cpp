#include "idresult_adv.h"

#define IDRESULT_EXIT_FAILURE 31

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
RESULT_SET* create_elem(const char *name, int age, const char *body) {
   RESULT_SET *ptr=(result_setp) malloc(sizeof(RESULT_SET));
   if (ptr!=NULL)
   {
   // allocates the memory for fields. Careful for chars, \0 must be counted
      ptr->name = (char*) malloc((strlen(name)+1)*sizeof(char));
      ptr->body = (char*) malloc((strlen(body)+1)*sizeof(char));
      if (ptr->name != NULL && ptr->body != NULL) 
      {
        // copy name and body 
        strcpy(ptr->name, name);
        strcpy(ptr->body, body);
        ptr->age=age;
        ptr->next=NULL;
      }
      else 
      {  // handle mem error
         free(ptr);
         ptr = NULL;
      }
   }
   return ptr;
}


//add on top
int add_elem(const char *name, int age, const char *body, RESULT_SET *head) {
   RESULT_SET *ptr=(result_setp) malloc(sizeof(RESULT_SET));
   if (ptr!=NULL)
   {
   // allocates the memory for fields. Careful for chars, \0 must be counted
      ptr->name = (char*) malloc((strlen(name)+1)*sizeof(char));
      ptr->body = (char*) malloc((strlen(body)+1)*sizeof(char));
      if (ptr->name != NULL && ptr->body != NULL) 
      {
        // copy name and body 
        strcpy(ptr->name, name);
        strcpy(ptr->body, body);
        ptr->age=age;
        ptr->next=head;
      }
      else 
      {  // handle mem error
         free(ptr);
         ptr = NULL;
      }
   }
   return (ptr != NULL);
}


//pop from top: display and remove
result_setp pop_elem(RESULT_SET *ptr) {
   if (ptr==NULL) {
     std::cout << "error (" << IDRESULT_EXIT_FAILURE << ")\n";
     return NULL;
   }
   else {
     result_setp pnext_c = ptr->next;
     printf("name: %s, age: %d, body: %s\n", ptr->name,ptr->age,ptr->body);
     free(ptr);
     return pnext_c;
   }
}

// remove particular element
int remove_elem(const char *matchname, RESULT_SET *ptr) {
   RESULT_SET* tmp = ptr;
   RESULT_SET* tmp_prev = NULL;
   int count=0;
   if (tmp==NULL) {
      std::cout << "empty linked list\n";
      return IDRESULT_EXIT_FAILURE;
   } else {
      while ((tmp != NULL) && (strcmp(tmp->name, matchname) != 0)) 
      {
         tmp_prev = tmp;
         tmp = tmp->next;
         count++;
      }
      if (tmp != NULL) {
         if (count == 0) {
            ptr = ptr->next;            
         } else {
            tmp_prev->next = tmp->next; //!!!
         }
         free(tmp);
      }
      else {
         std::cout << "name is not contained in this linked list\n";  
      }   
   }
   return 0;
}


//search for a particular cell for cell->name ="matchname"
// and replace its body and age
result_setp replace(char* matchname, int newage, char* newbody, int selected_id, RESULT_SET *head) {
   RESULT_SET* tmp;
   RESULT_SET* tmp_prev;
   tmp=head;
   tmp_prev=head;
   int count=0;
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
         RESULT_SET* tmp_ins = create_elem(tmp->name, newage, newbody);
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

//duplicate in reverse order
int copy_reverse_list(result_setp copylist, result_setp list) 
{   
   RESULT_SET* tmp=list;
   RESULT_SET* tmp_copy;
   int count=0;
   if (tmp==NULL) {
      printf("empty linked list");
      return -1;
   } else {
      while (tmp != NULL) {
         if (count==0) {
            tmp_copy = create_elem(tmp->name, tmp->age, tmp->body);
            tmp=tmp->next;
            count++;
         } else {
            add_elem(tmp->name, tmp->age, tmp->body, tmp_copy);
            tmp=tmp->next;
            count++;
         }       
      }
   }
   copylist = tmp_copy;
   return count;       
}  
   
/*
//search for a particular cell for cell->name ="matchname"
//and return another linked list with all matching cell for field cell->name
result_setp searchname(char* matchname, RESULT_SET *head) {
   RESULT_SET* tmp;
   RESULT_SET* tmp_prev;
   RESULT_SET* tmp_match;
   tmp=head;
   tmp_prev=head;
   int count=0;
   if (tmp==NULL) {
      printf("empty linked list");
      exit(1);
   } else {
      while ((tmp != NULL)) {
         if (strcmp(tmp->name, (const char*) matchname)) 
         {
            tmp_prev = tmp;
            tmp = tmp->next;
         }
         else { // if matching item
            if (count==0) { // matched item for the first time
 	       tmp_match = create(tmp->name, tmp->age, tmp->body);
               tmp_match->id =tmp->id;
               count++;
            }
            else { // successive matching items
               tmp_match = add(tmp->name, tmp->age, tmp->body, tmp_match);
	       tmp_match->id =tmp->id;
               count++;
	    }
         }//all in all the next record has to be examinated
         tmp_prev = tmp;
         tmp = tmp->next;
      }

   }
   return tmp_match;      
}
*/

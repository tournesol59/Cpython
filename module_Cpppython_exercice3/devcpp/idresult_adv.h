#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <iostream>

#ifdef __cplusplus
extern "C" 
{
#endif
    
// une petite structure de pile elementaire qui aide
// a resoudre la probleme de tableau de pointeur de pointeurs
typedef struct result_set RESULT_SET;
struct result_set {
    int id;
    char *name;
    int age;
    char *body;
    RESULT_SET *next;
};
typedef RESULT_SET* result_setp;

//methods

int list_length(result_setp head); 
 
RESULT_SET* create_elem(const char *name, int age, const char *body); 

int add_elem(const char *name, int age, const char *body, RESULT_SET *head);

result_setp pop_elem(RESULT_SET *ptr); 

int remove_elem(const char *matchname, RESULT_SET *ptr);

int copy_reverse_list(result_setp copylist, result_setp list);
/*
result_setp replace(char* matchname, int newage, char* newbody, int selected_id, RESULT_SET *head); 

result_setp searchname(char* matchname, RESULT_SET *head);
*/
#ifdef __cplusplus
}
#endif


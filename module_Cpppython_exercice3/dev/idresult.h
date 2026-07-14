#include <stdlib.h>
#include <string.h>
#include <stdio.h>

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

//methods
int list_length(result_setp head); 
 
RESULT_SET* create(char *name, int age, char *body); 
 
result_setp add(char *name, int age, char *body, RESULT_SET *head);
 
result_setp pop(RESULT_SET *ptr); 
 
result_setp replace(char* matchname, int newage, char* newbody, int selected_id, RESULT_SET *head); 
 
result_setp searchname(char* matchname, RESULT_SET *head);
 

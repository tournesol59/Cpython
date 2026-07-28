#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <iostream>

#ifdef __cplusplus
extern "C" 
{
#endif

// on garde la meme convention de nom que idresult_adv.h .cpp ce sont les noms des algorithmes
// et leur contenu (queue) qui changent
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

int queue_length(result_setp head);

RESULT_SET* queue_create(const char *name, int age, const char *body) ;

int queue_addlast(const char *name, int age, const char *body, RESULT_SET* head);

int queue_updatematch(const char *matchname, int newage, char* newbody, RESULT_SET *head);

#ifdef __cplusplus
}
#endif

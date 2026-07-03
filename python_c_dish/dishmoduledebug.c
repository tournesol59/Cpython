#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main (int argc, char **argv) {
    char command[60];
    int i,k,l,m;
    const char name[3]="sed";
    const char esc[2]=" ";
    //const char *format;
    //const char *input;
    //const char *flag;
    char flag[2];
    char format[10];
    char input[20];
    int sts;

    /*if (!PyArg_ParseTuple(args, "sss", &flag, &format, &input))
        return NULL;
    replaced by: */
    printf("entrez une option pour la commande sed (ex -i)");
    scanf("%s", flag);

    printf("entrez le pattern pour la command sed (ex s/H/h/g");
    scanf("%s", format);

    printf("entrez la chaine sur laquelle opere la commande sed (Hello..)");
    scanf("%s", input);

    k=strlen(flag);
    l=strlen(format);
    m=strlen(input);
    // concatenate in series into the "command" char
    // eg it will be similar to "sed  's/kong/hong' hong-kong
    strcpy(command, "");
    strcat(command, (char*) name);
    strcat(command, (char*) esc);
    strcat(command, (char*) flag);
    strcat(command, (char*) esc);
    strcat(command, (char*) format);
    strcat(command, (char*) esc);
    strcat(command, (char*) input);
    for (i=0; i < (60-3-1-k-1-l-1-m); i++) {
        strcat(command, (char*) esc);
    }
    // now we call shell to execute the command in string command:
    sts = system(command);
    return sts;  
}

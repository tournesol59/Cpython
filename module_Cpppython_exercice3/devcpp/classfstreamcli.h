#include <iostream>
#include <cstring>
#include <string>
//#include "idresultext.h"

//extern "C" int openforread(char *filename);
//extern "C" int handleregex(char *input, char *target);

class OPF_PRGMER {

    public:
       OPF_PRGMER(char* strconn, char* modus);
       OPF_PRGMER(const OPF_PRGMER &copyclass);
//       OPF_PRGMER &operator = (const OPF_PRGMER &copyclass);
       ~OPF_PRGMER();

       int openforreadlines();
       int writelineatend(char **charobject);
       
    private:
       char strconn[100];
       char modus[2];
       int count;
       int fserr;      

       int reportErrorResult(int errcode, char *errtext);
}; 

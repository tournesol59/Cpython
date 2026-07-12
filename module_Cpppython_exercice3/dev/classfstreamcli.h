#include <iostream>
#include <cstring>
//#include "idresultext.h"

extern "C" int openforread(char *filename);
extern "C" int handleregex(char *input, char *target);

class OPF_PRGMER {

    public:
       OPF_PRGMER(char* strconn, char* modus);
       OPF_PRGMER(const OPF_PRGMER &copyclass);
       OPF_PRGMER &operator = (const OPF_PRGMER &=copyclass);
       ~OPF_PRGMER();

    private:
       char strconn[100];
       char modus[2];
       int count;      
}

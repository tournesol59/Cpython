#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <sstream>
#include <assert.h>
#include "idresult_adv.h"

//extern "C" int openforread(char *filename);
//extern "C" int handleregex(char *input, char *target);

#define OPF_EXIT_OPENREAD 1
#define OPF_EXIT_OPENWRITE 2

class OPF_PRGMER {

    public:
       OPF_PRGMER(char* strconn, char* modus);
       OPF_PRGMER(const OPF_PRGMER &copyclass);
//       OPF_PRGMER &operator = (const OPF_PRGMER &copyclass);
       ~OPF_PRGMER();

       int openforreadlines();
       int openforreadlines_st(int nlines);
       int openforreadcharac();
       
       int writeresultsetToStr(std::string strhead, result_setp ptr, int nline);
       int writesetatend(result_setp resultLL_ptr);
       //int writelineatend(std::string strobj);
       //int writecharacatend(char **charobject);
       int readToResultSetHead(std::string strhead, result_setp resultLL_ptr, int iline);
       int readToResultSetData(std::string strline, result_setp resultLL_ptr, int iline);
       
    private:
       char strconn[100];
       char modus[2];
       int count;
       int fserr;      
       RESULT_SET* resultLL_p;
       
       int reportErrorResult(int errcode, std::string errtext);
}; 

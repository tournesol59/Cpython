// item definition and enum
#include <iostream>
#include <fstream>
#include "idjoinset.h"
#include <string>
#include <cstring>
#include <vector>
#include <regex.h>

using namespace std;

class OFT_PRGMER // open files/libelles of train transport
{
    public:
       OFT_PRGMER(const char* currentfilew, const char* modusw);
       OFT_PRGMER(const OFT_PRGMER &src);
       OFT_PRGMER & operator = (const OFT_PRGMER &src);
       ~OFT_PRGMER();
       
       int openforreadlines_st(int nlines);
       void interfaceUser();
       
       int findAll();
       int findByEmail(const char* email);
       int findByDepart(const char* depart);
       int findByDest(const char* destination);
       
       int ReportError(int errorToReport);
       int displayFirst();
       
    private:
       string m_currentfile;
       char m_modus[3];  
      // JoinSet_Ptr m_currentJSet; 
       vector<JOIN_SET> m_currentSets;
       // actually points to a single record of type JOIN_SET as first step 
       
};


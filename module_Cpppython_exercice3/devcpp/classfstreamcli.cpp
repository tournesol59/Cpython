#include <iostream>
#include <fstream>
//#include "idresultext.h"
#include "classfstreamcli.h"
#include <assert.h> 

/* declared in .h
class OPF_PRGMER {

    public:
       OPF_PRGMER(char* strconn, char* modus);
       OPF_PRGMER(const OPF_PRGMER &copyclass);
       OPF_PRGMER &operator = (const OPF_PRGMER &copyclass);
       ~OPF_PRGMER();

    private:
       char strconn[100];
       char modus[2];
       int count;      
}
       */

// constructor
OPF_PRGMER:: OPF_PRGMER(char *strconnect, char* modusw) {
   strcpy(strconn, (const char*) strconnect);
   strncpy(modus, (const char*) modusw, 1);
   modus[1]=0;
}

OPF_PRGMER:: OPF_PRGMER(const OPF_PRGMER &copyclass) {
   strcpy(strconn, (const char*) copyclass.strconn);
   strcpy(modus, (const char*) copyclass.modus);	
}
/*
OPF_PRGMER:: &operator= OPF_PRGMER (const OPF_PRGMER &copyclass) {
   strcpy(strconn, (const char*) copyclass.strconn);
   strcpy(modus, (const char*) copyclass.modus);	
}  	
*/
//destructor
OPF_PRGMER:: ~OPF_PRGMER() {}

//method opening a file and display its content
int OPF_PRGMER::openforreadlines() {
   std::string connfile(strconn);
   std::fstream fs;
   std::string line;

   fs.open(connfile, std::fstream::in);
   
   if (!fs.is_open())
   { 
      std::string errtext="could not open file";
      reportErrorResult(1, (char*) (errtext.c_str()));
/*
      fs.open(connfile, std::ios::out);
      fs << "Hello World" << std:endl;
      fs.close();
      fs.open(connfile); 
*/
   }
   while (getline(fs, line)) {
      std::cout << line << std::endl; 
   }
   fs.close();
   return 0; 
}

// write at the end 
int OPF_PRGMER::writelineatend(char **charobject) {
   int is_ok = 0;
   char **ptr=charobject;
   std::string line(*(ptr));
   int i=0;
   
   while ((ptr+i)!=NULL) {
      std::string obj(*(ptr+i++));
      line = line+" "+obj;  
   }
   assert(i==3); // debug dev code
   // now we have a line we can write to
   std::string fileoutname="sampoutput";
   std::fstream fs(fileoutname, std::ios::out);
   if (!fs.is_open()) {
      std::string errtext="cannnot open sampoutput file for writing";
      reportErrorResult(2, (char*) errtext.c_str());
   }
   else {
      fs << line << std::endl;
   }
   fs.close();
   is_ok=0;
   return is_ok;  
}



// ReportErrorResult
int OPF_PRGMER::reportErrorResult(int errcode, char *errtext) {
   std::string line;
   std::string msgtext(errtext);
  // std::ofstream fslog("Logfile", "w");
   
   // if (fslog) {
      switch(errcode) {
      	 case 1:
      		//fslog << "Fatal (1): file or directory not found : " << msgtext << std::endl;
      		std::cout << "Fatal (1): file or directory not found : " << msgtext << std::endl;
      		break; 
      	 case 2:
      		std::cout << "Fatal (2): fstream not enabled to create a new file : " << msgtext << std::endl;
      		break;
      	 }
      	 
   //}
  // fslog.close();
   return 0;
}




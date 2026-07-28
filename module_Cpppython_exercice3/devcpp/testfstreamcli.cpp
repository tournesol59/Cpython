#include <iostream>
#include "classfstreamcli.h"
//#include "idresult_adv.h"


//#define __cplusplus

int main() {
   char filesample[40];
   char modusw[2];
   int nline=3;
   // char** msgtext=(char**) malloc(2*sizeof (char*)); // must be followed by two other malloc of char*
   char msgtext[][40]={"Hello World", "Bonjour", ""}; // this has the advantage of setting a fixed length for each string

   result_setp resultLL_p;
   
   //FIRST TEST OF INIT A LINKEDLIST
   const char* name1="Sophie Bell"; 
   const char* body1="Medicine";
   resultLL_p=create_elem(name1, 1920, body1);
   const char* name2="Frances Allen";
   const char* body2="Mathematician";
   add_elem(name2, 1920, body2, resultLL_p);
   resultLL_p=pop_elem(resultLL_p);
   
   // READ CSV FILE
   strcpy(filesample, "sampledb");
   strcpy(modusw, "r");
   // read
   OPF_PRGMER opf_pgexample(filesample, modusw);
   opf_pgexample.openforreadlines_st(nline);
   
   //test RESULT_SET from C
   /*
      // Appel aux sous programmes
      
   */
      
   //write A TESTER
//   strcpy(msgtext[0], "This text file to print\n in text mode");
//   std::string msgexample(msgtext[0]);
//opf_pgexample.writelineatend(msgexample);
   
   
   //READ the file fileBin.bin file in binary mode
   //and WRITE the read buffer back in binary mode to another file
   //opf_example.writelineatend(msgtext);
   

   std::cout << "program testfstreamcli normally exited\n";
   return 0;
}

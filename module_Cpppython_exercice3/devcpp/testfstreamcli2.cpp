#include <iostream>
#include "classfstreamcli.h"
//#include "idresult_adv.h"


//#define __cplusplus

int main(int argc, char** argv) {
   char filesample[40];
   char modusw[2];
   int nlines=atol(argv[1]);  // 3;
   // char** msgtext=(char**) malloc(2*sizeof (char*)); // must be followed by two other malloc of char*
   //char msgtext[][40]={"Hello World", "Bonjour", ""}; // this has the advantage of setting a fixed length for each string

   result_setp resultLL_p;
   
   //FIRST TEST OF INIT A LINKEDLIST: it is either "elem" (meaning LIFO) or "queue" (FIFO)
   const char* name1="Sophie Bell"; 
   const char* body1="Medicine";
   resultLL_p=create_elem(name1, 1920, body1);
   const char* name2="Frances Allen";
   const char* body2="Mathematician";
   add_elem(name2, 1920, body2, resultLL_p);
   resultLL_p=pop_elem(resultLL_p);
   
   // READ CSV FILE with dedicated class, classfstreamcli
   strcpy(filesample, "sampledb");
   strcpy(modusw, "r");
   // read
   OPF_PRGMER opf_pgexample(filesample, modusw);
   opf_pgexample.openforreadlines_st(nlines);
   
   //test RESULT_SET from C
   /*
      // Appel aux sous programmes
      
   */
      
   //write A TESTER
    opf_pgexample.writesetatend(resultLL_p);

   std::cout << "program testfstreamcli normally exited\n";
   return 0;
}

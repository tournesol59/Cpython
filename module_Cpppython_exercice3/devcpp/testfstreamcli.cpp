#include <iostream>
#include "classfstreamcli.h"

int main() {
   char filesample[40];
   char modusw[2];
   char** msgtext=(char**) malloc(2*sizeof (char*));
   
   strcpy(filesample, "sampledb");
   strcpy(modusw, "r");
   
   OPF_PRGMER opf_pgex(filesample, modusw);
   opf_pgex.openforreadlines();
   
   strcpy(*msgtext, "Hello World");
   opf_pgex.writelineatend(msgtext);
   
   std::cout << "program testfstreamcli normally exited\n";
   return 0;
}

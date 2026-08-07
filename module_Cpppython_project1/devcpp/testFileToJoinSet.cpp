#include <iostream>
#include "classfstreamcli.h"
//#include "idresult_adv.h"


//#define __cplusplus

int main(int argc, char** argv) {
   char filesample[40];
   char modusw[2];
   int nlines=atol(argv[1]);  // 3;

   // READ CSV FILE with dedicated class, classfstreamcli
   strcpy(filesample, "sampledba");
   strcpy(modusw, "r");
   
   // read
   OFT_PRGMER open_prog(filesample, modusw);
   open_prog.openforreadlines_st(nlines);
   
   //here we call the function display of a single struct JOIN_SET
   return 0;
}

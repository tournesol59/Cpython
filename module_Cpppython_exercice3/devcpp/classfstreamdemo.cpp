// Copy a file
#include <iostream>
#include <fstream>
#include "classfstreamcli.h"
/*
int OPF_PRGMER::writecharacatend(char** charobject) {
   // pour un premier test: recopie
   std::ifstream infile("fileBin.bin", std::ifstream::binary);
   std::ofstream outfile("fileOut.bin", std::ofstream::binary);

   // get size of file
   infile.seekg(0, infile.end);
   long size = infile.tellg();
   infile.seekg(0);

   // allocate memory for file content
   char* buffer = new char[size];

   // read content of infile
   infile.read(buffer, size);

   // write to outfile
   outfile.write(buffer, size);

   // release dynamically-allocated memory
   delete[] buffer;

   outfile.close();
   infile.close();
   return 0;
}
*/
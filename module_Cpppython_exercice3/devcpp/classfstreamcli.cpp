#include <iostream>
#include <fstream>
//#include "idresult.h"
#include "classfstreamcli.h"
#include <assert.h> 
#include <regex>

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
   resultLL_p=create_elem("",20,"");
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
   // open in text
   fs.open(connfile, std::fstream::in);
   
   if (!fs.is_open())
   { 
      std::string errtext="could not open file";
      reportErrorResult(OPF_EXIT_OPENREAD, errtext);

   }
   // read a line after another until EOF
   while (getline(fs, line)) {
      std::cout << line << std::endl; 
   }
   fs.close();
   return 0; 
}

// another method, this time parse each line assumed in csv format and save in member list
int OPF_PRGMER::openforreadlines_st(int nlines) {
   int is_ok=0;
   std::string connfile(strconn);
   std::ifstream fs;
   std::string line;
   std::string strhead("name,age,function,");
   int count=0;
   char name[40];
   int age;
   char body[40];
   
   // setting the streams exception mask, before calling open()
   fs.exceptions(std::ifstream::failbit);
   try {
      fs.open(connfile, std::fstream::in);
       // inside the try, double check
      if (!fs.is_open())
      {
         std::string errtext="could not open file";
         reportErrorResult(1, errtext);
      }
   // traverse the file, from which we assume to be formed of 
   // comma separators : fixed 3 commas (one comma at the eol also)

      // create a structure of linked list
      RESULT_SET* element_list=create_elem("", 20, ""); // arbitrary values for first
   
      while (getline(fs, line) && (count<nlines)) {  // adding a counter was necessary
         // because strangely the handling of the std::ifstream::failbit causes error 
         // when the last line is reached
      //while (getline(fs, line)) {
          
         // use a stringstream object to split each line with the separator (delim)
         std::istringstream iss(line);
         std::string temp;
         char delim=',';
         std::string x;
         int i=0;
         
         if (count == 0) {
            readToResultSetHead(strhead, element_list, count);
            count++;
         }
   	  	 else {
         // else
         // readToResultSetData(line, element_list, count); // this proc shall be created
         	std::cout << "Parsing line " << std::endl << line << std::endl;
         	const std::regex e1("[a-zA-Z]+");
         	const std::regex e2("[0-9]+");
         	while (getline(iss, temp, delim)) {
            	if (std::istringstream(temp) >> x) {
               	switch(i)
               	{
                    case 0:
                        if (std::regex_match(x, e1)) {  
                    	    strcpy(name, (const char*) x.c_str());
                    	}
                    	else {
                    	   memset(name, '\0', sizeof name); 
                    	}
                    	break;   
                  	case 1:
                  	    //if (std::regex_match(x, e2)) {
                  	       
                    	age = atol(x.c_str());
                    	break;       
                  	case 2:
                        if (std::regex_match(x, e1)) {  
                    	   strcpy(body, (const char*) x.c_str());
                    	}
                    	else {
                    	   memset(body, '\0', sizeof body); 
                    	}
                    	break;                                      
               	}
               	i++;
            }   
         }
         assert(i==3);
         if ((name[0]=='\0') || (body[0]=='\0')) {
            std::cout << "formatting of the fields of line " << count << " aren t good \n";   
         }
         
         // add a new cell at the top of the linked list with last infos
         add_elem((const char*) name, age, (const char*) body, element_list);
         i=0; // for next line 
         count++;
         }//end if count==0 (wrong identation)
      }
   
      // This last step create a local copy as a class member variable
//copy_reverse_list(resultLL_p, element_list);
      // this has an advantage: pop order in list resultLL_p correspond to order in the file 
      //check
      fs.close();
      return is_ok;
      
   }  catch (const std::exception& e) {
      std::ostringstream msg;
      msg << "Opening file '" << connfile << "' failed, it either does not exist or is not accessible.";
      reportErrorResult(1, msg.str());
      throw std::runtime_error(msg.str()); 
      is_ok=1;
      return is_ok;  
   }
   
}

int OPF_PRGMER::openforreadcharac() {
   int is_ok = 0;        
   std::string str;  
   // open the binary file    
   std::ifstream file("fileBin.bin", std::ios::binary);
   if (!file) {
      std::cerr << "Error opening or reading\n";
      return 1;
   }
   // Read the length of the string(size)
   size_t strLength;
   file.read(reinterpret_cast<char*>(&strLength), sizeof(strLength));
   
   // Allocate memory for the string and read the data
   char* buffer = new char[strLength+1];
   file.read(buffer, strLength);
   // Null-terminate the string 
   buffer[strLength] = '\0';
   // Convert the buffer
   str = buffer;
   // Print file data
   std::cout << "File Data: " << str;
   delete[] buffer;
   
   // Close file
   file.close();
   return is_ok;   
}


int OPF_PRGMER::readToResultSetHead(std::string strhead, result_setp resultLL_ptr, int iline)
{
   // strhead contains the titles of the columns, only commas-separated strings
   // stritem contains a line of different values
   // resultLL_ptr contains a pointer to the linked list output structure
   
   int is_ok=1;
   //faire avec strtok
   std::string sub(","); // keep it for counting inside a std::string ..
   char sep[]=",";       // .. but parsing will be used with C the strtok function
   char strcopy[100];
   //char crep='\0'; //
     // copy the header input string into a char[]
   strcpy(strcopy, strhead.c_str());
      
   size_t count=0; // for better and safe behavior
// count the occurences of separator sub
   for (size_t offset=strhead.find(sub); offset != std::string::npos; offset=strhead.find(sub, offset+sub.length()))
   {
      ++count;
   }
   std::cout << "number of fields separated by a comma: " <<  count << std::endl;

   //finally we decided that the better option to parse the string separated by commas
   // would be the C call to strtok, indeed it can be modified for other separator or multiple separators
   
   std::vector<std::string> strhead_split;
   char *token;
   std::string stritem;
   
   size_t idx = 0;
   token = strtok(strcopy, sep);
   while ((token != NULL) && (idx < count)) {  //bug fixe mais pas entierement
      // while there are token in string
      printf(" %s\n", token);
      stritem.assign(token);
      strhead_split.push_back(stritem);
      
      token = strtok(NULL, sep);
      idx++;   
   }
   assert(idx==3);
   for (std::vector<std::string>::iterator it=strhead_split.begin(); it != strhead_split.end(); ++it) {
       std::cout << *it << "\t";   
   }
   std::cout << std::endl;
   
   
   // next we should parse a line of content: this time we cannot use a vector because the type of 
   // items have different types, in the same line
   is_ok=0;
   return is_ok;
}

// this proc shall be used for parsing a single line but of Data (not header)
int readToResultSetData(std::string strline, result_setp resultLL_ptr, int iline) {
   int is_ok=0;
   //...
   return is_ok;
}
          
int OPF_PRGMER::writeresultsetToStr(std::string strhead, result_setp ptr, int nline) {
   // concatenate a cell object data into a single string
   
   int is_ok=1;
   std::string sub(","); // keep it for counting inside a std::string ..
  // char sep[]=",";       // .. but parsing will be used with C the strtok function
  // char strcopy[100];
   //char crep='\0'; // to delete
   //std::string rep(crep); // to delete
 
   
   /* to delete
   // now operate the decumposition of the header string into a vector of sub-strings   
   std::string sub(",");
   
   std::vector<std::string> strhead_split;
   std::string stritem;
   for (size_t offset=strhead.find(sub); offset != std::string::npos; offset=strhead.find(sub, offset+sub.length()))
   {
      stritem = strhead.substr(offset, offset+sub.length()-1);
      std::cout << stritem;
      strhead_split.push_back(stritem);
   }
   */
   is_ok=0;
   return is_ok;
}

int OPF_PRGMER::writesetatend(result_setp resultLL_ptr)  {
   int is_ok=1;
   int iline=0;
   std::string strheader = "name,age,function";
   //result_setp resultLL_ptr = queue_create("Frances Allen", 1920, "Mathematician");
   
   // test spliter
   is_ok=writeresultsetToStr(strheader, resultLL_ptr, iline);
   
   return is_ok;
}


/*
// internal subroutine to be used by writelineatend for demo puposes of string treatment
std::string printTextReplaced(std::string str, std::string rep) {
   //sub= "\n";
   std::string sub("text");
   //std::string rep= "example";

   std::cout << "original: " << str << "\n"; // origin
   int count=0; // for debug
   for (size_t offset=str.find(sub); offset != std::string::npos; offset=str.find(sub, offset+rep.length()))
   {
      ++count;
   }
   
   for (size_t offset=str.find(sub); offset != std::string::npos; offset=str.find(sub, offset+rep.length()))
   {
      str.replace(offset, sub.length(), rep);   
   }
   
   std::cout << "replaced " << count << " words in: " << str << "\n"; // end printout

   return str;
}

// write at the end 
//int OPF_PRGMER::writelineatend(char *charobject[12]) { // ne plus jamais faire ca
int OPF_PRGMER::writelineatend(std::string strobj) // preferer les string
{
   int is_ok = 0;
   
   /// a bannir
   //while ((ptr+i)!=NULL) {
   //   std::string obj(*(ptr+i++));
   //   line = line+" "+obj;  
   //}
   //
   // count eol
   int count;
   std::string subeol("\n");
   for (size_t offset=strobj.find(subeol); offset != std::string::npos;
         offset = strobj.find(subeol, offset+subeol.length()))
   {
     ++count;
   }

   //int i=0;
   // size_t offset=strobj.find(subeol);
   //.. a finir, on regarde tout d'abord le replacement simple
   
   std::string line=printTextReplaced(strobj, "example");
      
      
   // now we have a line we can write to
   std::string fileoutname="sampoutput";
   std::fstream fs(fileoutname, std::ios::out);
   if (!fs.is_open()) {
      std::string errtext="cannot open sampoutput file for writing";
      reportErrorResult(OPF_EXIT_OPENWRITE, errtext);
   }
   else {
      fs << line << std::endl;
   }
   fs.close();
   is_ok=0;
   return is_ok;  
}
*/


// ReportErrorResult
int OPF_PRGMER::reportErrorResult(int errcode, std::string errtext) {
   std::string line;
   std::string msgtext(errtext);
  // std::ofstream fslog("Logfile", "w");
   
   // if (fslog) {
      switch(errcode)
      {
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




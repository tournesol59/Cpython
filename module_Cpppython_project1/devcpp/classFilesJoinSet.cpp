#include "classFilesJoinSet.h"

using namespace std;

OFT_PRGMER::OFT_PRGMER(const char* currentfilew, const char* modusw) :
   m_currentfile(currentfilew)
   {
   	 strcpy(m_modus, modusw);     
   }
	
OFT_PRGMER::OFT_PRGMER(const OFT_PRGMER &src) :
   m_currentfile(src.m_currentfile)
   {
     strcpy(m_modus, (const char*) src.m_modus);	
   }
	
//OFT_PRGMER::&operator=OFT_PGRMER(const OFT_PRGMER &src) :
//   m_currentfile(src.m_currentfile)
//   {
//      strcpy(m_modus, (const char*) src.m_modus);	
//   }
   
   //destructor
OFT_PRGMER:: ~OFT_PRGMER() {}

//method, parse each line assumed in csv format and save in member list
int OFT_PRGMER::openforreadlines_st(int nlines) 
{
   int is_ok=0;
   std::ifstream fs;
   std::string line;
   std::string strhead("email,insurance,refnum,departure,destination,spare,spareref");
   int count=0;
   char email[40];
   int insuran_nr;
   char refnum[40];
   char departure[40];
   char destination[40];
   char spare[40];
   int spareref;
   
   // setting the streams exception mask, before calling open()
   fs.exceptions(std::ifstream::failbit);
   try {
      fs.open(m_currentfile, std::fstream::in);
       // inside the try, double check
       // inside the try, double check
      if (!fs.is_open())
      {
         std::string errtext="could not open file";
         // TBC preview a reportGenerator
      }
      
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
            //readToResultSetHead(strhead, line, element_list, count); // 
            cout << strhead << endl;
            count++;
         }
   	  	 else {         // else
         // readToResultSetData(line, element_list, count); // this proc shall be created
         	std::cout << "Parsing line " << std::endl << line << std::endl;
         //	const std::regex e1("[a-zA-Z]+"); // erreur ici: does not name a type
         //	const std::regex e2("[0-9]+");
         	while (getline(iss, temp, delim)) {
            	if (std::istringstream(temp) >> x) {
            	   switch(i)
               	   {
                    case 0:
                       strcpy(email, x.c_str());
                       break;
                    case 1:
                       insuran_nr = (int) ((Insurance) atol(x));
                       break;
                    case 2:
                       strcpy(refnum, x.c_str());
                       break;
                    case 3:
                       strcpy(departure, x.c_str());
                       break;
                    case 4:
                       strcpy(destination, x).c_str();
                       break;
                    case 5:
                       strcpy(spare, x.c_str());
                       break;
                    case 6:
                       spareref = (int) atol(x);
                       break; 
            	   }
            	}	
            }
            m_currentJet = create_joinset(const char* email, int insurance, const char* refnum, const char* departure, const char* destination, const char* string spare, int spareref);
            count++;
      }
      
      fs.close();
      return is_ok;
      
   }  catch (const std::exception& e) {
      std::ostringstream msg;
      msg << "err 1 Opening file '" << m_currentfile << "' failed";
      ReportError(1);
      throw std::runtime_error(msg.str()); 
      is_ok=1;
      return is_ok;  
   }
}
int OFT_PRGMER::displayFirst() {
   // TBC
   vector<JOIN_SET>::iterator it=m_currentSets.begin();
   
      return 0;   
}


/*
void OFT_PRGMER::interfaceUser()
{
   string x1;
   int choice_request_1;
   //int choice_request_2;
   string x2a;
   string x2b;
   string choice_field2;
   string choice_value2;
   vector<string> field_names = {
      "email", "insurance", "refnum", "departure", "destination", "spare1", "spare2"
   }
   
   cout << "Select a type of request\n";
   cin >> x1;
   choice_request_1 = atol(x1);
   switch (choice_request_1) {
      case 0:
         if (findAll()) return 1;
         break;
      case 1:
         cout << "Select a field name\n";
         cin >> x2a;
         choice_field2 = x2a;
         cout << "Select a field value\n";
         cin >> x2b;
         choice_value2 = x2b;
         // validation
         int count_field2=0;
         for (vector<string>::iterator it=field_names.begin() ; it != field_names.end(); ; ++it) 
         {
            if (choice_field2 == *it) break;
            count_field2++ 
         }
         if (count_field < 7) {
               switch(count_field2) {
                  case 0:
                     if (findByEmail((const char*) choice_value2)) return 2; 
                     break; 
                  case 3:
                     if (findByDepart((const char*) choice_value2)) return 2; 
                     break; 
                  case 4:
                     if (findByDest((const char*) choice_value2)) return 2; 
                     break; 
               }
          } 
          break;
   }
}
  */ 
int OFT_PRGMER::ReportError(int errorToReport) {
   LookupErrorMessage errorMessage = LookupErrorMessage(errorToReport);
   std::cout << errorMessage.getText() << std::endl;
  
   return 0; 
}

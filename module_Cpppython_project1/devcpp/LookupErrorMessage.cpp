#include "LookupErrorMessage"
using namespace std;

LookupErrorMessage::LookupErrorMessage(int errorToReport) :
   m_code({1}),
   m_messages({"temp"})
   {
      errorMessageInit();   
      // insert code to store the message appropriate to m_code==errorToReport
   }

void LookupErrorMessage::ErrorMessageInit() 
{
   int code;
   string message;   
   // replace the following by a list of messages
   code = ERR_CODE_1;
   message("error opening stream");
   m_codes.push_back(code);
   m_messages.push_back(message);
   
   code = ERR_CODE_2;
   message("error writing stream");
   m_codes.push_back(code);
   m_messages.push_back(message);
}

   


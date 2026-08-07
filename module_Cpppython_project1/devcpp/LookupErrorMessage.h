#include <iostream>
#include <vector>

#define ERR_CODE_1 1
#define ERR_CODE_2 2

using namespace std;

class LookupErrorMessage {
   public:
      void ErrorMessageInit();
      LookupErrorMessage(int errorToReport);
      ~LookupErrorMessage();
      string getText();

   private:
      m_errorToReport;
      m_text;
      vector<int> m_codes;
      vector<string> m_messages;
};
